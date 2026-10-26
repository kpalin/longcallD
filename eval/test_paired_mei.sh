#!/usr/bin/env bash
#
# Paired-mode mobile element insertions. A somatic insertion carrying a TSD and
# a polyA tail may be called on --som-mei-alt tumour reads, as -s mode has
# always allowed, provided the normal holds no alt read at all. The generic
# indel gates -- min_alt_dp, the indel VAF floor and the per-base artifact
# checks -- used to leave every such call Ambiguous.
#
# The data are synthetic (eval/make_paired_mei.py): error-free reads with exact
# CIGARs, so each case has exactly the alt reads it says it has.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
TE=$ROOT/anno/AluY_L1_SVA_cons_noPA.fa
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

pass=0; fail=0
chk() {
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

python3 -I "$ROOT/eval/make_paired_mei.py" "$ROOT/test_data/chr11_2M.fa" "$TE" "$T" || exit 1
samtools faidx "$T/ref.fa"
for s in tumor normal; do
  samtools sort -o "$T/$s.bam" "$T/$s.sam" 2>/dev/null && samtools index "$T/$s.bam"
done

"$BIN" call -t 1 -T "$TE" "$T/ref.fa" "$T/tumor.bam" -N "$T/normal.bam" > "$T/out.vcf" 2> "$T/log" \
  || { echo "longcallD failed:"; tail "$T/log"; exit 1; }

# The insertion record for event $1: SV INS within 50 bp of its site.
rec() {
  local p len; read -r p len < <(awk -v e="$1" '$1==e {print $2, $3}' "$T/events.tsv")
  awk -F'\t' -v p="$p" -v n="$len" '!/^#/ && $2>p-50 && $2<p+50 && length($5)-length($4)>n-10' "$T/out.vcf" | head -1
}
has() { [[ "$1" == *"$2"* ]] && echo yes || echo no; }

a=$(rec A); b=$(rec B); c=$(rec C); d=$(rec D)
chk "A (2 tumour reads, clean normal) is reported"   yes "$([[ -n $a ]] && echo yes || echo no)"
chk "A is SOMATIC"                                   yes "$(has "$a" 'SOMATIC;')"
chk "A FILTER is PASS"                               PASS "$(cut -f7 <<< "$a")"
chk "A is annotated MEI +AluY"                       yes "$(has "$a" 'REPNAME=+AluY')"
chk "B (1 tumour read) is SOMATIC"                   yes "$(has "$b" 'SOMATIC;')"
chk "B is annotated MEI -AluY"                       yes "$(has "$b" 'REPNAME=-AluY')"
chk "C is reported (so the next check is not vacuous)" yes "$([[ -n $c ]] && echo yes || echo no)"
chk "C (1 alt read in the normal) is not SOMATIC"    no  "$(has "$c" 'SOMATIC;')"
chk "D is reported"                                  yes "$([[ -n $d ]] && echo yes || echo no)"
chk "D (no TSD/polyA) is not SOMATIC"                no  "$(has "$d" 'SOMATIC;')"
chk "D carries no MEI annotation"                    no  "$(has "$d" 'MEI;')"

# --som-mei-alt raises the bar back above B's single read.
"$BIN" call -t 1 -T "$TE" --som-mei-alt 2 "$T/ref.fa" "$T/tumor.bam" -N "$T/normal.bam" > "$T/out2.vcf" 2>/dev/null
b2=$(awk -F'\t' -v p="$(awk '$1=="B"{print $2}' "$T/events.tsv")" '!/^#/ && $2>p-50 && $2<p+50 && length($5)>300' "$T/out2.vcf" | head -1)
chk "--som-mei-alt 2: B is no longer SOMATIC"        no  "$(has "$b2" 'SOMATIC;')"
a2=$(awk -F'\t' -v p="$(awk '$1=="A"{print $2}' "$T/events.tsv")" '!/^#/ && $2>p-50 && $2<p+50 && length($5)>300' "$T/out2.vcf" | head -1)
chk "--som-mei-alt 2: A is still SOMATIC"            yes "$(has "$a2" 'SOMATIC;')"

echo
if (( fail == 0 )); then echo "   paired MEI: all $pass assertions passed"; exit 0
else echo "   paired MEI: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
