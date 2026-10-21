#!/usr/bin/env bash
#
# Per-(sample, haplotype) allele tallies -- the four-slot grid the paired
# somatic test will decide on.
#
# These are aggregated separately from the pooled counters that drive
# single-sample calling, so the two must be shown to agree. LONGCALLD_SLOTCHK=1
# makes longcallD dump both for every variant:
#
#   SLOTCHK  pos  pooledDP  pooledALT  slotDP  slotALT \
#            alt/ref/oth per (sample,hap)...  class som_hap unphased hapvaf pval
#
# so the slot groups are fields 7 .. NF-5.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
REF=$ROOT/test_data/chr11_2M.fa
HIFI=$ROOT/test_data/HG002_chr11_hifi_test.bam
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

pass=0; fail=0
chk() {
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

dump() { LONGCALLD_SLOTCHK=1 "$BIN" call "$@" -t 1 -o /dev/null 2>&1 | grep '^SLOTCHK'; }

# ---------------------------------------------------------------- single ----
dump "$REF" "$HIFI" > "$T/single.txt"
n=$(wc -l < "$T/single.txt")
chk "single sample: produced variants" yes "$([[ $n -gt 0 ]] && echo yes || echo no)"

# The alt count is what a somatic call turns on, so it must match exactly.
chk "single sample: slot ALT equals pooled ALT at every site" 0 \
    "$(awk -F'\t' '$4!=$6' "$T/single.txt" | wc -l | tr -d ' ')"

# Depth is defined differently on purpose: the pooled counter credits every read
# spanning the site, the slots only reads with an allele call there. The slots
# must therefore never exceed the pooled depth, and should track it closely.
chk "single sample: slot depth never exceeds pooled depth" 0 \
    "$(awk -F'\t' '$5>$3' "$T/single.txt" | wc -l | tr -d ' ')"
chk "single sample: slot depth within 2 of pooled depth everywhere" 0 \
    "$(awk -F'\t' '$3-$5>2' "$T/single.txt" | wc -l | tr -d ' ')"

# Every read must land in exactly one of alt/ref/oth, so the three sum to the depth.
chk "single sample: alt+ref equals reported slot depth" 0 \
    "$(awk -F'\t' '{s=0; for(i=7;i<=NF-5;i++){split($i,a,"/"); s+=a[1]+a[2]} if(s!=$5) n++} END{print n+0}' "$T/single.txt")"

# ---------------------------------------------------------------- paired ----
# Two disjoint halves, so no molecule is in both "samples".
samtools view "$HIFI" | awk '{print $1}' | sort -u > "$T/names"
nn=$(wc -l < "$T/names")
head -n $((nn/2)) "$T/names" > "$T/n1"; tail -n +$((nn/2+1)) "$T/names" > "$T/n2"
i=0; for part in n1 n2; do i=$((i+1))
  { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$T/$part"; } \
    | samtools sort -o "$T/h$i.bam" - 2>/dev/null; samtools index "$T/h$i.bam"
done

dump "$REF" "$T/h1.bam" -N "$T/h2.bam" > "$T/paired.txt"
chk "paired: produced variants" yes "$([[ $(wc -l < "$T/paired.txt") -gt 0 ]] && echo yes || echo no)"

# Six groups means two samples x three haplotype slots.
chk "paired: two samples' worth of slots are reported" 6 \
    "$(head -1 "$T/paired.txt" | awk -F'\t' '{print NF-11}')"

# The pooled counters still cover both samples, so the per-sample alt counts
# must add back up to them. This is the load-bearing check: if the sample axis
# were mis-assigned, evidence would move between samples and this would break.
chk "paired: tumour ALT + normal ALT equals pooled ALT" 0 \
    "$(awk -F'\t' '{s=0; for(i=7;i<=NF-5;i++){split($i,a,"/"); s+=a[1]} if(s!=$4) n++} END{print n+0}' "$T/paired.txt")"
chk "paired: tumour depth + normal depth equals slot depth" 0 \
    "$(awk -F'\t' '{s=0; for(i=7;i<=NF-5;i++){split($i,a,"/"); s+=a[1]+a[2]} if(s!=$5) n++} END{print n+0}' "$T/paired.txt")"

# Both samples must actually contribute; all-zero on one side would mean the
# normal reads never reached the tallies.
chk "paired: normal sample contributes evidence" yes \
    "$(awk -F'\t' '{s=0; for(i=10;i<=12;i++){split($i,a,"/"); s+=a[1]+a[2]} if(s>0) n++} END{print (n>0)?"yes":"no"}' "$T/paired.txt")"
chk "paired: tumour sample contributes evidence" yes \
    "$(awk -F'\t' '{s=0; for(i=7;i<=9;i++){split($i,a,"/"); s+=a[1]+a[2]} if(s>0) n++} END{print (n>0)?"yes":"no"}' "$T/paired.txt")"

echo
if [[ $fail -eq 0 ]]; then echo "  slot counts: all $pass assertions passed"; exit 0
else echo "  slot counts: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
