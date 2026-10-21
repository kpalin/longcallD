#!/usr/bin/env bash
#
# E2/E3: score longcallD against the SEQC2 HCC1395 somatic truth set on chr1,
# the same configuration the DeepSomatic paper reports for PacBio HiFi.
#
# Two modes:
#   -m subtract   call each sample separately and subtract the normal's calls
#                 from the tumour's. This is the baseline the paper uses for
#                 DeepVariant, and needs no paired-calling code -- it is the
#                 "before" number the paired caller has to beat.
#   -m paired     longcallD paired (default)
#
# Published chr1 PacBio HiFi figures for the same truth set and region, for
# comparison (Park et al. 2025, Nat Biotechnol, Supplementary Tables 2-3):
#
#       caller                      SNV F1    indel F1
#       DeepSomatic single-cancer   0.9536    0.8151
#       ClairS v0.2.0               0.9310    0.4550
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
DATA=${DATA:-/home/ubuntu/castle_tn}
OUT=${OUT:-$DATA/bench}
THREADS=${THREADS:-4}
MODE=paired
REGION=${REGION:-chr1}
EXTRA=${EXTRA:-}

TUMOR=$DATA/HCC1395_pacbio.tumor.chr1.bam
NORMAL=$DATA/HCC1395_pacbio.normal.chr1.bam
REF=$DATA/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna

usage() {
  cat <<EOF
Usage: $0 [-m paired|subtract] [-r REGION] [-t THREADS] [-o OUTDIR] [-x "EXTRA ARGS"]

  -m  calling mode (default: $MODE)
  -r  region to call; the scoring is always restricted to chr1 (default: $REGION)
  -t  threads (default: $THREADS)
  -o  output directory (default: $OUT)
  -x  extra arguments passed through to longcallD

Needs the fixtures from eval/fetch_data.sh and the SIF from eval/build_happy_sif.sh.
EOF
  exit 1
}
while getopts "m:r:t:o:x:h" o; do
  case $o in
    m) MODE=$OPTARG ;; r) REGION=$OPTARG ;; t) THREADS=$OPTARG ;;
    o) OUT=$OPTARG ;; x) EXTRA=$OPTARG ;; *) usage ;;
  esac
done

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

for f in "$BIN" "$TUMOR" "$NORMAL" "$REF"; do
  [[ -s $f ]] || die "missing $f (run eval/fetch_data.sh?)"
done
mkdir -p "$OUT"

case $MODE in
  paired)
    log "longcallD paired on $REGION ($THREADS threads)"
    /usr/bin/time -v "$BIN" paired -t "$THREADS" $EXTRA "$REF" "$TUMOR" "$NORMAL" "$REGION" \
        -o "$OUT/paired.vcf" > "$OUT/paired.stdout" 2> "$OUT/paired.log" \
      || die "longcallD paired failed; see $OUT/paired.log"
    # som.py scores PASS records only, so the somatic set is what PASSes and
    # carries the SOMATIC flag. Germline calls also PASS and must come out.
    bcftools view -f PASS -i 'INFO/SOMATIC=1' -Oz -o "$OUT/query.vcf.gz" "$OUT/paired.vcf" 2>/dev/null \
      || die "could not extract somatic calls"
    ;;
  subtract)
    log "single-sample calls, then subtraction (the DeepVariant-style baseline)"
    for s in tumor normal; do
      bam=$([[ $s == tumor ]] && echo "$TUMOR" || echo "$NORMAL")
      "$BIN" call -s -t "$THREADS" $EXTRA "$REF" "$bam" "$REGION" -o "$OUT/$s.vcf" 2> "$OUT/$s.log" \
        || die "longcallD call failed on $s; see $OUT/$s.log"
      bgzip -f "$OUT/$s.vcf"; tabix -f -p vcf "$OUT/$s.vcf.gz"
    done
    bcftools isec -C -w1 -Oz -o "$OUT/query.vcf.gz" "$OUT/tumor.vcf.gz" "$OUT/normal.vcf.gz" 2>/dev/null \
      || die "bcftools isec failed"
    ;;
  *) usage ;;
esac

tabix -f -p vcf "$OUT/query.vcf.gz"
log "$(bcftools view -H "$OUT/query.vcf.gz" | wc -l) candidate somatic calls"

log "scoring with som.py"
"$ROOT/eval/sompy.sh" -q "$OUT/query.vcf.gz" -o "$OUT/sompy.$MODE" -l chr1 \
  || die "som.py failed"

cat <<EOF

Reference figures for the same truth set, region and technology
(Park et al. 2025, Nat Biotechnol):

  DeepSomatic single-cancer   SNV F1 0.9536   indel F1 0.8151
  ClairS v0.2.0               SNV F1 0.9310   indel F1 0.4550

EOF
log "done: $OUT/sompy.$MODE.stats.csv"
