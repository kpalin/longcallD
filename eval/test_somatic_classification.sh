#!/usr/bin/env bash
#
# The four-slot somatic test: a variant is somatic only if its alt occupies
# exactly one of (N,H1) (N,H2) (T,H1) (T,H2), and that slot is a tumour one.
#
# The bundled files are a single sample, so splitting them in half gives a
# tumour/normal pair that shares a germline and contains no somatic variants --
# a negative control, and the strongest specificity statement available here.
# Positive controls are manufactured with eval/inject_somatic.py, which writes a
# new base into one haplotype of one half (partitioning reads by their allele at
# a nearby het site) and rewrites the =/X CIGAR to match.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
REF=$ROOT/test_data/chr11_2M.fa
HIFI=$ROOT/test_data/HG002_chr11_hifi_test.bam
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

# A well-covered het SNV to tell the two haplotypes apart, a nearby
# variant-free site to inject into, and the base to write there.
HET=1258256; HET_ALT=G; SITE=1258300; NEW=A
# 1e-3 cannot be reached with the handful of alt reads this shallow test file
# provides -- the smallest attainable p for 4 alt reads here is 1.3e-3 -- so the
# threshold is relaxed for the test. Real data at 50x+ does not need this.
PVAL=0.01

pass=0; fail=0
chk() {
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

samtools view "$HIFI" | awk '{print $1}' | sort -u > "$T/names"
n=$(wc -l < "$T/names")
head -n $((n/2)) "$T/names" > "$T/n1"; tail -n +$((n/2+1)) "$T/names" > "$T/n2"

# $1 = read-name list, $2 = output bam, $3 = injection fraction ("" = none)
build() {
  if [[ -z ${3:-} ]]; then
    { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$1"; } \
      | samtools sort -o "$2" - 2>/dev/null
  else
    { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$1"; } \
      | python3 -I "$ROOT/eval/inject_somatic.py" "$HET" "$HET_ALT" "$SITE" "$NEW" "$3" 2>/dev/null \
      | samtools sort -o "$2" - 2>/dev/null
  fi
  samtools index "$2"
}
build "$T/n1" "$T/tumor_clean.bam"
build "$T/n2" "$T/normal_clean.bam"
build "$T/n1" "$T/tumor_clonal.bam" 1.0
build "$T/n1" "$T/tumor_sub.bam" 0.5
build "$T/n2" "$T/normal_inj.bam" 1.0

slots() { LONGCALLD_SLOTCHK=1 "$BIN" call --som-pval "$PVAL" "$REF" "$1" -N "$2" -t 1 -o /dev/null 2>&1 | grep '^SLOTCHK'; }
class_at() { awk -F'\t' -v p="$2" '$2==p {print $13}' "$1"; }
field_at() { awk -F'\t' -v p="$2" -v c="$3" '$2==p {print $c}' "$1"; }
count_class() { awk -F'\t' -v c="$2" '$13==c' "$1" | wc -l | tr -d ' '; }

SOMATIC=1; GERMLINE=2; NORMAL_ONLY=3

# ------------------------------------------------- negative control --------
slots "$T/tumor_clean.bam" "$T/normal_clean.bam" > "$T/neg.txt"
chk "negative control: variants were classified" yes \
    "$([[ $(wc -l < "$T/neg.txt") -gt 0 ]] && echo yes || echo no)"
chk "negative control: no somatic call anywhere" 0 "$(count_class "$T/neg.txt" $SOMATIC)"

# ------------------------------------------------- clonal somatic ----------
slots "$T/tumor_clonal.bam" "$T/normal_clean.bam" > "$T/clonal.txt"
chk "clonal: the injected site is called somatic" "$SOMATIC" "$(class_at "$T/clonal.txt" $SITE)"
chk "clonal: it is the only somatic call" 1 "$(count_class "$T/clonal.txt" $SOMATIC)"
chk "clonal: attributed to a single haplotype" yes \
    "$([[ "$(field_at "$T/clonal.txt" $SITE 14)" =~ ^[12]$ ]] && echo yes || echo no)"
chk "clonal: decided on the slot grid, not the pooled fallback" 0 "$(field_at "$T/clonal.txt" $SITE 15)"
# A clonal variant occupies its whole haplotype, so the within-haplotype VAF is
# 1.0 -- where the pooled VAF would read about 0.25 and look subclonal.
chk "clonal: HAPVAF is 1.0" "1.0000" "$(field_at "$T/clonal.txt" $SITE 16)"
chk "clonal: everything else is classified as before" \
    "$(awk -F'\t' -v p=$SITE '$2!=p {print $13}' "$T/clonal.txt" | sort | uniq -c | md5sum)" \
    "$(awk -F'\t' '{print $13}' "$T/neg.txt" | sort | uniq -c | md5sum)"

# ------------------------------------------------- subclonal somatic -------
slots "$T/tumor_sub.bam" "$T/normal_clean.bam" > "$T/sub.txt"
chk "subclonal: still called somatic" "$SOMATIC" "$(class_at "$T/sub.txt" $SITE)"
# Injected into half the haplotype's reads, so HAPVAF must drop well below 1.0
# while staying clearly above zero -- it tracks cellular fraction.
chk "subclonal: HAPVAF tracks the injected fraction" yes \
    "$(awk -v v="$(field_at "$T/sub.txt" $SITE 16)" 'BEGIN{print (v>0.3 && v<0.8) ? "yes" : "no"}')"

# ------------------------------------------------- normal-only -------------
# The user-facing requirement: a variant present only in the germline is
# suspect, and must never be reported as somatic.
slots "$T/tumor_clean.bam" "$T/normal_inj.bam" > "$T/normonly.txt"
chk "normal-only: flagged as normal-only" "$NORMAL_ONLY" "$(class_at "$T/normonly.txt" $SITE)"
chk "normal-only: not called somatic" 0 "$(awk -F'\t' -v p=$SITE '$2==p && $13==1' "$T/normonly.txt" | wc -l | tr -d ' ')"

echo
if [[ $fail -eq 0 ]]; then echo "  somatic classification: all $pass assertions passed"; exit 0
else echo "  somatic classification: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
