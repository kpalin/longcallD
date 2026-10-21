#!/usr/bin/env bash
#
# E1 gate -- prove eval/sompy.sh is wired correctly before any longcallD number
# is believed.
#
# DeepSomatic does not publish the query VCF behind its case-study table (the
# brain-genomics-public and deepvariant buckets carry truth sets, BAMs and
# panels of normals only), so we cannot re-derive their 3324/147/249 directly.
# Instead we check the two things that table actually certifies:
#
#   A. total.truth must be 3573 on chr1, split 3440 SNV / 133 indel -- exactly
#      the published column. This pins the truth VCF, the HC BED, -l chr1, the
#      reference and bcftools norm to the paper's configuration.
#   B. FP/FN accounting, by perturbing the truth by a known amount and
#      requiring som.py to report precisely that amount.
#
set -euo pipefail

DATA=${DATA:-/home/ubuntu/castle_tn}
TRUTH=${TRUTH:-$DATA/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}
OUTDIR=${OUTDIR:-$DATA/eval_out/validate}
HERE=$(cd "$(dirname "$0")" && pwd)

# Published chr1 truth composition (DeepSomatic PacBio case study).
EXP_TOTAL=3573 EXP_SNV=3440 EXP_INDEL=133
N_DROP=100 N_MUTATE=50

log()  { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die()  { log "FAIL: $*"; exit 1; }
mkdir -p "$OUTDIR"

# field from a som.py stats.csv: $1=file $2=type $3=column
stat_of() { python3 -c "
import csv,sys
for r in csv.DictReader(open(sys.argv[1])):
    if r['type']==sys.argv[2]:
        print(int(float(r[sys.argv[3]]))); break
" "$1" "$2" "$3"; }

check() { # name expected actual
  if [[ "$2" == "$3" ]]; then log "  ok   $1 = $3"
  else die "$1: expected $2, got $3"; fi
}

CHR1="$OUTDIR/truth.chr1.vcf.gz"
if [[ ! -s $CHR1 ]]; then
  log "subsetting truth to chr1"
  bcftools view -r chr1 -Oz -o "$CHR1" "$TRUTH"
  tabix -f -p vcf "$CHR1"
fi

# ---------------------------------------------------------------- test A ----
log "A: identity (truth vs itself)"
"$HERE/sompy.sh" -q "$CHR1" -o "$OUTDIR/a_identity" >/dev/null 2>&1 \
  || die "som.py failed on test A"
S="$OUTDIR/a_identity.stats.csv"
check "A total.truth"      "$EXP_TOTAL" "$(stat_of "$S" records total.truth)"
check "A SNV total.truth"  "$EXP_SNV"   "$(stat_of "$S" SNVs    total.truth)"
check "A indel total.truth" "$EXP_INDEL" "$(stat_of "$S" indels total.truth)"
check "A tp" "$EXP_TOTAL" "$(stat_of "$S" records tp)"
check "A fp" 0            "$(stat_of "$S" records fp)"
check "A fn" 0            "$(stat_of "$S" records fn)"

# ---------------------------------------------------------------- test B ----
# Drop N_DROP SNVs from the query. Pure false negatives, no false positives.
log "B: drop $N_DROP SNVs -> expect fn=$N_DROP, fp=0"
Q="$OUTDIR/b_dropped.vcf"
bcftools view -h "$CHR1" > "$Q"
bcftools view -H "$CHR1" | awk -v n="$N_DROP" '
  BEGIN{d=0}
  { if (length($4)==1 && length($5)==1 && d<n) { d++; next } print }
' >> "$Q"
"$HERE/sompy.sh" -q "$Q" -o "$OUTDIR/b_dropped" >/dev/null 2>&1 \
  || die "som.py failed on test B"
S="$OUTDIR/b_dropped.stats.csv"
check "B fn" "$N_DROP"                       "$(stat_of "$S" records fn)"
check "B fp" 0                               "$(stat_of "$S" records fp)"
check "B tp" "$((EXP_TOTAL - N_DROP))"       "$(stat_of "$S" records tp)"

# ---------------------------------------------------------------- test C ----
# Change the ALT base on N_MUTATE SNVs. POS/REF stay valid so the records
# survive `bcftools norm -c x` and the HC BED, but the ALT no longer matches --
# som.py should score each as one FP and one FN, since matching is exact
# CHROM/POS/REF/ALT string equality.
log "C: mutate ALT on $N_MUTATE SNVs -> expect fp=$N_MUTATE, fn=$N_MUTATE"
Q="$OUTDIR/c_mutated.vcf"
bcftools view -h "$CHR1" > "$Q"
bcftools view -H "$CHR1" | awk -v n="$N_MUTATE" 'BEGIN{OFS="\t"; m=0}
  {
    if (length($4)==1 && length($5)==1 && m<n) {
      split("A C G T", b, " ")
      for (i=1;i<=4;i++) if (b[i]!=$4 && b[i]!=$5) { $5=b[i]; break }
      m++
    }
    print
  }' >> "$Q"
"$HERE/sompy.sh" -q "$Q" -o "$OUTDIR/c_mutated" >/dev/null 2>&1 \
  || die "som.py failed on test C"
S="$OUTDIR/c_mutated.stats.csv"
check "C fp" "$N_MUTATE"                   "$(stat_of "$S" records fp)"
check "C fn" "$N_MUTATE"                   "$(stat_of "$S" records fn)"
check "C tp" "$((EXP_TOTAL - N_MUTATE))"   "$(stat_of "$S" records tp)"

log "PASS: som.py harness validated (identity, false negatives, false positives)"
