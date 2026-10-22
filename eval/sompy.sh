#!/usr/bin/env bash
#
# E1 -- score a somatic VCF against the SEQC2 HCC1395 high-confidence truth set
# using the EXACT som.py the DeepSomatic paper used (pkrusche/hap.py v0.3.9,
# 2017), run from an apptainer SIF.
#
# som.py reports recall and precision but not F1; the paper's F1 column is
# derived downstream, so this script computes it from the .stats.csv.
#
# Note on obtaining the SIF: `apptainer pull docker://pkrusche/hap.py` FAILS.
# Every tag of that image is a 2017-era Docker manifest schema v1, which
# apptainer's go-containerregistry backend cannot read. skopeo still reads
# schema v1, so convert first (see eval/build_happy_sif.sh).
#
set -euo pipefail

DATA=${DATA:-/home/ubuntu/castle_tn}
SIF=${SIF:-$DATA/happy.sif}
TRUTH=${TRUTH:-$DATA/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}
HCBED=${HCBED:-$DATA/High-Confidence_Regions_v1.2.bed}
REF=${REF:-$DATA/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna}
LOC=${LOC:-chr1}
TARGETS=${TARGETS:-}

usage() {
  cat <<EOF
Usage: $0 -q QUERY.vcf[.gz] -o OUT_PREFIX [-l LOCATION] [-T TARGETS.bed]

  -q   query (caller) VCF. Only PASS records are scored by som.py.
  -o   output prefix; som.py writes <prefix>.stats.csv and <prefix>.features.csv
  -l   location restriction (default: $LOC). Applied IN ADDITION to the HC BED.
  -T   capture-target BED, for exome evaluation (default: none)

Env: DATA SIF TRUTH HCBED REF
EOF
  exit 1
}

QUERY="" OUT=""
while getopts "q:o:l:T:h" o; do
  case $o in
    q) QUERY=$OPTARG ;;
    o) OUT=$OPTARG ;;
    l) LOC=$OPTARG ;;
    T) TARGETS=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -n $QUERY && -n $OUT ]] || usage

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

for f in "$SIF" "$TRUTH" "$HCBED" "$REF" "$QUERY"; do
  [[ -s $f ]] || die "missing or empty: $f"
done
command -v apptainer >/dev/null || die "apptainer not found"

# Restrict the query to small variants before scoring.
#
# som.py preprocesses with a shell pipeline that ends
#   bcftools view -R <HC BED> ... | bcftools norm -f ref -c x -D
# and `view -R` emits a record ONCE PER OVERLAPPING BED INTERVAL. A deletion
# long enough to span several high-confidence intervals therefore comes out
# several times, at a position earlier than the intervals that follow, so the
# stream runs backwards and som.py dies indexing its own intermediate:
#
#   [E::hts_idx_push] unsorted positions on sequence #1: 194483236 followed by
#   194482054
#
# The culprit on chr1 is a 3,129 bp deletion in a TA microsatellite, emitted
# four times. This is a limitation of the 2017 tooling, not something the
# caller can avoid, and it cannot be worked around by normalising first --
# the duplication happens inside som.py, after our file is read.
#
# Excluding these records is also the correct comparison. The SEQC2 small
# variant truth set tops out at 55 bp, so an SV-sized call can never be a true
# positive; and the published figures we compare against come from small
# variant callers, which do not emit such records at all. The bound is applied
# to every arm, and the count dropped is reported so its effect stays visible.
MAXLEN=${MAXLEN:-60}
SMALL="${OUT}.small.vcf.gz"
bcftools view -i "strlen(REF)<=$MAXLEN && strlen(ALT)<=$MAXLEN" -Oz -o "$SMALL" "$QUERY" 2>/dev/null \
  || die "failed to restrict $QUERY to small variants"
tabix -f -p vcf "$SMALL"
n_in=$(bcftools view -H "$QUERY" | wc -l); n_out=$(bcftools view -H "$SMALL" | wc -l)
log "small-variant bound (<= ${MAXLEN} bp): kept $n_out of $n_in records, dropped $((n_in - n_out))"
QUERY="$SMALL"

SCRATCH="${OUT}.scratch"; rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"

# Bind every distinct parent directory the inputs live in. apptainer maps $HOME
# and the cwd automatically, but data held elsewhere is invisible without this.
binds=()
while read -r d; do binds+=(--bind "$d:$d"); done < <(
  for f in "$TRUTH" "$HCBED" "$REF" "$QUERY" "$SCRATCH" "${TARGETS:-$REF}"; do
    dirname "$(readlink -f "$f")"
  done | sort -u
)

args=(-N "$TRUTH" "$QUERY" -r "$REF" -o "$OUT"
      --feature-table generic -R "$HCBED"
      --scratch-prefix "$SCRATCH")
[[ -n $LOC ]] && args+=(-l "$LOC")
[[ -n $TARGETS ]] && args+=(-T "$TARGETS")

log "som.py ${args[*]}"
apptainer exec "${binds[@]}" "$SIF" /opt/hap.py/bin/som.py "${args[@]}"

[[ -s ${OUT}.stats.csv ]] || die "som.py produced no ${OUT}.stats.csv"

# Add the F1 column the paper reports but som.py does not emit.
python3 - "$OUT.stats.csv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
want = ("indels", "SNVs", "records")
hdr = f"{'type':<9} {'total.truth':>11} {'total.query':>11} {'tp':>6} {'fp':>6} {'fn':>6} {'recall':>8} {'precision':>9} {'F1':>8}"
print(hdr); print("-" * len(hdr))
for r in rows:
    if r.get("type") not in want:
        continue
    tp, fp, fn = (int(r[k]) for k in ("tp", "fp", "fn"))
    rec = tp / (tp + fn) if tp + fn else 0.0
    pre = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * rec * pre / (rec + pre) if rec + pre else 0.0
    print(f"{r['type']:<9} {int(r['total.truth']):>11} {int(r['total.query']):>11} "
          f"{tp:>6} {fp:>6} {fn:>6} {rec:>8.4f} {pre:>9.4f} {f1:>8.4f}")
PY

rm -rf "$SCRATCH"
log "done: ${OUT}.stats.csv ${OUT}.features.csv"
