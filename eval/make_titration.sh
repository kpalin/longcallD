#!/usr/bin/env bash
#
# E4: tumour-purity titration, built the way the DeepSomatic paper builds it.
#
# A tumour arm at purity p is p of the tumour's reads plus (1-p) of the
# normal's, at a constant total depth so recall is not confounded with coverage.
#
# The normal is split into two disjoint halves first: one is spiked into the
# tumour arm as contamination, the other is the matched normal. No molecule is
# ever both the contamination and the evidence against it -- without that split
# the same read can appear in both arms and the contamination is invisible to
# the caller.
#
# Depths are measured, not assumed: the achievable total is bounded by the
# spike-in half, since the lowest-purity point needs the most normal reads.
# A requested fraction above 1.0 is a hard error rather than a silently
# shallower arm, which would quietly turn a purity curve into a depth curve.
#
set -uo pipefail

DATA=${DATA:-/home/ubuntu/castle_tn}
OUT=${OUT:-$DATA/titration}
REGION=${REGION:-chr1}
THREADS=${THREADS:-4}
PURITIES=${PURITIES:-"1.00 0.90 0.75 0.50 0.25"}
TUMOR=${TUMOR:-$DATA/HCC1395_pacbio.tumor.chr1.bam}
NORMAL=${NORMAL:-$DATA/HCC1395_pacbio.normal.chr1.bam}
HCBED=${HCBED:-$DATA/High-Confidence_Regions_v1.2.bed}
# samtools view -s takes SEED.FRACTION in one number. Distinct seeds per draw,
# as the paper's scripts do, so the subsets are independent.
SEED_T=11; SEED_N=42

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
for f in "$TUMOR" "$NORMAL" "$HCBED"; do [[ -s $f ]] || die "missing $f"; done
mkdir -p "$OUT"

# Mean depth over the high-confidence regions of one contig, sampled every
# 50th position -- plenty for a depth estimate and far cheaper than all of it.
depth_of() {
  local bam=$1 bed=$OUT/.hc.$REGION.bed
  [[ -s $bed ]] || awk -v c="${REGION%%:*}" '$1==c' "$HCBED" > "$bed"
  samtools depth -b "$bed" -Q 20 "$bam" 2>/dev/null \
    | awk 'NR%50==0 {t+=$3; n++} END{if(n) printf "%.2f", t/n; else print 0}'
}

# Split the normal into two disjoint halves by read name.
if [[ ! -s $OUT/normal_eval.bam || ! -s $OUT/normal_spike.bam ]]; then
  log "splitting the normal into disjoint evaluation and spike-in halves"
  samtools view "$NORMAL" | awk '{print $1}' | sort -u > "$OUT/.names"
  n=$(wc -l < "$OUT/.names")
  head -n $((n/2)) "$OUT/.names" > "$OUT/.n_eval"; tail -n +$((n/2+1)) "$OUT/.names" > "$OUT/.n_spike"
  for part in eval spike; do
    { samtools view --no-PG -H "$NORMAL"; samtools view "$NORMAL" | grep -Fwf "$OUT/.n_$part"; } \
      | samtools sort -@2 -o "$OUT/normal_$part.bam" - 2>/dev/null
    samtools index "$OUT/normal_$part.bam"
  done
  shared=$(comm -12 <(samtools view "$OUT/normal_eval.bam" | cut -f1 | sort -u) \
                    <(samtools view "$OUT/normal_spike.bam" | cut -f1 | sort -u) | wc -l)
  [[ $shared -eq 0 ]] || die "the two normal halves share $shared read names"
fi

d_tumor=$(depth_of "$TUMOR")
d_spike=$(depth_of "$OUT/normal_spike.bam")
d_eval=$(depth_of "$OUT/normal_eval.bam")
log "measured depth over HC regions: tumour ${d_tumor}x, spike-in half ${d_spike}x, matched normal ${d_eval}x"

# The binding constraint is the lowest purity: it needs (1-p_min) x TOTAL from
# the spike-in half. Leave 5% headroom for the depth estimate itself.
p_min=$(echo "$PURITIES" | tr ' ' '\n' | sort -g | head -1)
TOTAL=${TOTAL:-$(python3 -c "print(min(float('$d_tumor'), float('$d_spike')/(1-$p_min+1e-9))*0.95)")}
log "tumour-arm total depth held at $(printf '%.1f' "$TOTAL")x across the series"

for p in $PURITIES; do
  out=$OUT/tumor_p${p}.bam
  if [[ -s $out ]] && samtools quickcheck -q "$out" 2>/dev/null; then log "purity $p: already built"; continue; fi
  ft=$(python3 -c "f=$p*$TOTAL/float('$d_tumor'); print('%.4f'%f)")
  fn=$(python3 -c "f=(1-$p)*$TOTAL/float('$d_spike'); print('%.4f'%f)")
  python3 -c "
import sys
ft, fn = float('$ft'), float('$fn')
if ft > 1.0: sys.exit('tumour fraction %.3f > 1 at purity $p: not enough tumour depth' % ft)
if fn > 1.0: sys.exit('spike-in fraction %.3f > 1 at purity $p: not enough normal depth' % fn)
" || die "purity $p is unreachable from the available depth"

  log "purity $p: ${ft} of the tumour + ${fn} of the spike-in half"
  samtools view -@2 -b -s "${SEED_T}.${ft#0.}" "$TUMOR" > "$OUT/.t.bam" 2>/dev/null
  if [[ $(python3 -c "print(1 if $fn > 0 else 0)") == 1 ]]; then
    samtools view -@2 -b -s "${SEED_N}.${fn#0.}" "$OUT/normal_spike.bam" > "$OUT/.n.bam" 2>/dev/null
    samtools merge -@2 -f "$out" "$OUT/.t.bam" "$OUT/.n.bam" 2>/dev/null || die "merge failed at purity $p"
    rm -f "${OUT:?}/.n.bam"
  else
    mv "$OUT/.t.bam" "$out"
  fi
  rm -f "${OUT:?}/.t.bam"
  samtools index "$out"
  log "purity $p: $(samtools view -c "$out") records, $(depth_of "$out")x"
done

log "done. tumour arms in $OUT, matched normal is $OUT/normal_eval.bam"
