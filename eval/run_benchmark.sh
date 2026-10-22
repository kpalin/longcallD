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
# -t 2, not 4. longcallD batches n_threads*4 regions of 500 kb and keeps
# several batches in flight, so the thread count sets the size of the live read
# set. Measured peak RSS, paired, chr1 from position 1: 1.5 GB at 10 Mb, 2.5 GB
# at 25 Mb, 3.2 GB at 50 Mb -- growth with region length, sublinear but real.
# At -t 4 that doubles, and whole-chr1 reached 6.7 GB and was OOM-killed on this
# 7 GB box.
THREADS=${THREADS:-2}
# Call in windows rather than one pass, which bounds the live set and makes the
# run resumable: a window whose VCF already exists is skipped. 0 disables.
WINDOW=${WINDOW:-25000000}
# Hard ceiling per invocation, so a blow-up fails this run with an allocation
# error instead of inviting the OOM killer to choose a victim.
MEMCAP_KB=${MEMCAP_KB:-5000000}
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
  -w  window size in bp for chunked, resumable calling; 0 disables (default: $WINDOW)

Needs the fixtures from eval/fetch_data.sh and the SIF from eval/build_happy_sif.sh.
EOF
  exit 1
}
while getopts "m:r:t:o:x:w:h" o; do
  case $o in
    m) MODE=$OPTARG ;; r) REGION=$OPTARG ;; t) THREADS=$OPTARG ;; w) WINDOW=$OPTARG ;;
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
    call_window() { # $1 = region string, $2 = output vcf
      ( ulimit -v "$MEMCAP_KB"
        "$BIN" paired -t "$THREADS" $EXTRA "$REF" "$TUMOR" "$NORMAL" "$1" -o "$2" ) 2>> "$OUT/paired.log"
    }
    if [[ ${WINDOW:-0} -gt 0 && $REGION != *:* ]]; then
      ctg=$REGION
      len=$(awk -v c="$ctg" '$1==c {print $2}' "$REF.fai" 2>/dev/null)
      [[ -n $len ]] || die "no length for $ctg in $REF.fai"
      log "longcallD paired on $ctg in $(( (len + WINDOW - 1) / WINDOW )) windows of ${WINDOW} bp ($THREADS threads)"
      : > "$OUT/paired.log"; parts=(); failed_windows=()
      for ((beg=1; beg<=len; beg+=WINDOW)); do
        end=$(( beg + WINDOW - 1 )); (( end > len )) && end=$len
        part="$OUT/part.$(printf '%012d' "$beg")_${end}.vcf"; parts+=("$part")
        if [[ -s $part ]]; then log "  $ctg:$beg-$end already done"; continue; fi
        log "  $ctg:$beg-$end"
        if call_window "$ctg:$beg-$end" "$part.tmp"; then
          mv "$part.tmp" "$part"
        else
          # Keep going and say so at the end. One window that cannot be called
          # should cost its own variants, not the other 90% of the chromosome.
          log "  WARNING: $ctg:$beg-$end failed; continuing without it"
          failed_windows+=("$ctg:$beg-$end"); rm -f "$part.tmp"
          : > "$part"   # placeholder so a resumed run does not retry it blindly
        fi
      done
      if [[ ${#failed_windows[@]} -gt 0 ]]; then
        log "${#failed_windows[@]} window(s) produced no calls: ${failed_windows[*]}"
        log "their truth variants will count as false negatives, which is the honest accounting"
      fi
      log "concatenating ${#parts[@]} windows"
      hdr=""; for f in "${parts[@]}"; do [[ -s $f ]] && { hdr=$f; break; }; done
      [[ -n $hdr ]] || die "every window failed; nothing to score"
      { grep '^#' "$hdr"; for f in "${parts[@]}"; do [[ -s $f ]] && grep -v '^#' "$f"; done; } > "$OUT/paired.vcf"
    else
      log "longcallD paired on $REGION ($THREADS threads)"
      : > "$OUT/paired.log"
      call_window "$REGION" "$OUT/paired.vcf" || die "longcallD paired failed; see $OUT/paired.log"
    fi
    # som.py scores PASS records only, so the somatic set is what PASSes and
    # carries the SOMATIC flag. Germline calls also PASS and must come out.
    bcftools view -f PASS -i 'INFO/SOMATIC=1' -Oz -o "$OUT/query.vcf.gz" "$OUT/paired.vcf" 2>/dev/null \
      || die "could not extract somatic calls"
    ;;
  subtract)
    log "single-sample calls, then subtraction (the DeepVariant-style baseline)"
    for s in tumor normal; do
      bam=$([[ $s == tumor ]] && echo "$TUMOR" || echo "$NORMAL")
      : > "$OUT/$s.log"
      # Windowed for the same reason as the paired branch: memory grows with
      # region length, and single-sample whole-chr1 is large enough to matter
      # even though it is lighter than paired.
      if [[ ${WINDOW:-0} -gt 0 && $REGION != *:* ]]; then
        len=$(awk -v c="$REGION" '$1==c {print $2}' "$REF.fai")
        [[ -n $len ]] || die "no length for $REGION in $REF.fai"
        sparts=()
        for ((beg=1; beg<=len; beg+=WINDOW)); do
          end=$(( beg + WINDOW - 1 )); (( end > len )) && end=$len
          sp="$OUT/$s.part.$(printf '%012d' "$beg").vcf"; sparts+=("$sp")
          [[ -s $sp ]] && continue
          log "  $s $REGION:$beg-$end"
          ( ulimit -v "$MEMCAP_KB"
            "$BIN" call -s -t "$THREADS" $EXTRA "$REF" "$bam" "$REGION:$beg-$end" -o "$sp.tmp" ) 2>> "$OUT/$s.log" \
            && mv "$sp.tmp" "$sp" || { log "  WARNING: $s $REGION:$beg-$end failed; continuing"; rm -f "$sp.tmp"; : > "$sp"; }
        done
        hdr=""; for f in "${sparts[@]}"; do [[ -s $f ]] && { hdr=$f; break; }; done
        [[ -n $hdr ]] || die "every $s window failed"
        { grep '^#' "$hdr"; for f in "${sparts[@]}"; do [[ -s $f ]] && grep -v '^#' "$f"; done; } > "$OUT/$s.vcf"
      else
        ( ulimit -v "$MEMCAP_KB"; "$BIN" call -s -t "$THREADS" $EXTRA "$REF" "$bam" "$REGION" -o "$OUT/$s.vcf" ) 2>> "$OUT/$s.log" \
          || die "longcallD call failed on $s; see $OUT/$s.log"
      fi
      bgzip -f "$OUT/$s.vcf"; tabix -f -p vcf "$OUT/$s.vcf.gz"
    done
    bcftools isec -C -w1 -Oz -o "$OUT/query.vcf.gz" "$OUT/tumor.vcf.gz" "$OUT/normal.vcf.gz" 2>/dev/null \
      || die "bcftools isec failed"
    ;;
  *) usage ;;
esac

tabix -f -p vcf "$OUT/query.vcf.gz"
log "$(bcftools view -H "$OUT/query.vcf.gz" | wc -l) candidate somatic calls"

# The somatic F1 is the headline, but the other classes carry the specificity
# story: with a real matched normal there should be few normal-only calls, and
# the grid/pooled split says how much of the answer rests on the weaker test.
if [[ -s $OUT/paired.vcf ]]; then
  log "call composition:"
  awk -F'\t' '!/^#/ {
        c = ($8 ~ /(^|;)SOMATIC(;|$)/) ? "SOMATIC" :
            ($8 ~ /(^|;)LOH(;|$)/)     ? "LOH" :
            ($8 ~ /NORMALONLY/)        ? "NORMAL_ONLY" : "germline/other"
        n[c]++; if (c=="SOMATIC") { if ($8 ~ /UNPHASED/) pooled++; else grid++ }
      } END {
        for (k in n) printf "    %-15s %d\n", k, n[k]
        printf "    somatic decided on the haplotype grid: %d, on pooled counts: %d\n", grid+0, pooled+0
      }' "$OUT/paired.vcf" >&2
fi

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
