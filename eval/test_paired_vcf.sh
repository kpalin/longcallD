#!/usr/bin/env bash
#
# The two-sample VCF: both columns, the paired verdict in FILTER/INFO, and the
# evidence grid in FORMAT.
#
# The load-bearing claim it has to support is that phasing is shared between
# tumour and normal -- same phase sets, same orientation -- because the whole
# four-slot test rests on haplotype h meaning the same chromosome copy in both.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
REF=$ROOT/test_data/chr11_2M.fa
HIFI=$ROOT/test_data/HG002_chr11_hifi_test.bam
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
HET=1258256; HET_ALT=G; SITE=1258300; NEW=A; PVAL=0.01

pass=0; fail=0
chk() {
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

samtools view "$HIFI" | awk '{print $1}' | sort -u > "$T/names"
n=$(wc -l < "$T/names"); head -n $((n/2)) "$T/names" > "$T/n1"; tail -n +$((n/2+1)) "$T/names" > "$T/n2"
build() {
  if [[ -z ${3:-} ]]; then
    { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$1"; } | samtools sort -o "$2" - 2>/dev/null
  else
    { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$1"; } \
      | python3 -I "$ROOT/eval/inject_somatic.py" "$HET" "$HET_ALT" "$SITE" "$NEW" "$3" 2>/dev/null \
      | samtools sort -o "$2" - 2>/dev/null
  fi
  samtools index "$2"
}
build "$T/n1" "$T/tumor.bam" 1.0
build "$T/n2" "$T/normal.bam"
build "$T/n1" "$T/tumor_clean.bam"
build "$T/n2" "$T/normal_inj.bam" 1.0

"$BIN" call --som-pval $PVAL "$REF" "$T/tumor.bam" -N "$T/normal.bam" -t 1 -o "$T/paired.vcf" 2>/dev/null
"$BIN" call --som-pval $PVAL "$REF" "$T/tumor_clean.bam" -N "$T/normal_inj.bam" -t 1 -o "$T/normonly.vcf" 2>/dev/null

# ------------------------------------------------------------- structure ---
# Column names come from each sample's own @RG SM. Both halves of the bundled
# file carry SM:HG002, so they collide and are suffixed rather than emitted
# twice -- bcftools rejects a duplicate-sample VCF outright.
chk "two sample columns, normal first" "HG002_N HG002_T" \
    "$(grep '^#CHROM' "$T/paired.vcf" | cut -f10,11 | tr '\t' ' ')"
chk "the column names are distinct" distinct \
    "$(grep '^#CHROM' "$T/paired.vcf" | awk '{print ($10==$11) ? "same" : "distinct"}')"
chk "bcftools accepts the header" 0 "$(bcftools view -h "$T/paired.vcf" >/dev/null 2>&1; echo $?)"
for tag in NormalOnly LowDP Ambiguous; do
  chk "FILTER $tag is declared" 1 "$(grep -c "##FILTER=<ID=$tag," "$T/paired.vcf")"
done
# UNPHASED is INFO, not FILTER: as a filter it made every such call non-PASS, so
# som.py dropped it and the two evidence regimes could not be compared.
chk "UNPHASED is INFO, not FILTER" "0 1" \
    "$(grep -c '##FILTER=<ID=Unphased,' "$T/paired.vcf") $(grep -c '##INFO=<ID=UNPHASED,' "$T/paired.vcf")"
for tag in NORMALONLY LOH SOMHAP HAPVAF SPV; do
  chk "INFO $tag is declared" 1 "$(grep -c "##INFO=<ID=$tag," "$T/paired.vcf")"
done
for tag in HDP HAD; do
  chk "FORMAT $tag is declared" 1 "$(grep -c "##FORMAT=<ID=$tag," "$T/paired.vcf")"
done

# --------------------------------------------------------------- somatic ---
rec=$(grep -v '^#' "$T/paired.vcf" | awk -v p=$SITE '$2==p')
chk "the injected somatic variant is emitted" yes "$([[ -n $rec ]] && echo yes || echo no)"
chk "it passes the filter" "PASS" "$(echo "$rec" | cut -f7)"
chk "it is flagged SOMATIC" yes "$(echo "$rec" | cut -f8 | grep -q 'SOMATIC' && echo yes || echo no)"
chk "it names the haplotype" yes "$(echo "$rec" | cut -f8 | grep -qE 'SOMHAP=[12]' && echo yes || echo no)"
chk "it reports HAPVAF and SPV" yes \
    "$(echo "$rec" | cut -f8 | grep -q 'HAPVAF=' && echo "$rec" | cut -f8 | grep -q 'SPV=' && echo yes || echo no)"
# The normal must show no alt on either haplotype; the tumour must show it on one.
chk "normal has no alt on either haplotype" "0,0" \
    "$(bcftools query -f '%POS\t[%HAD\t]\n' "$T/paired.vcf" 2>/dev/null | awk -v p=$SITE '$1==p{print $2}')"
chk "tumour has alt on exactly one haplotype" yes \
    "$(bcftools query -f '%POS\t[%HAD\t]\n' "$T/paired.vcf" 2>/dev/null | awk -v p=$SITE '$1==p{split($3,a,","); print ((a[1]>0)!=(a[2]>0)) ? "yes" : "no"}')"

# ----------------------------------------------------------- normal-only ---
# The stated requirement: present only in the germline means flagged and
# excluded from PASS, but still inspectable rather than silently dropped.
rec2=$(grep -v '^#' "$T/normonly.vcf" | awk -v p=$SITE '$2==p')
chk "a normal-only variant is still emitted" yes "$([[ -n $rec2 ]] && echo yes || echo no)"
chk "it is filtered NormalOnly" "NormalOnly" "$(echo "$rec2" | cut -f7)"
chk "it carries the NORMALONLY flag" yes "$(echo "$rec2" | cut -f8 | grep -q 'NORMALONLY' && echo yes || echo no)"
chk "bcftools -f PASS drops it" 0 \
    "$(bcftools view -f PASS "$T/normonly.vcf" 2>/dev/null | grep -v '^#' | awk -v p=$SITE '$2==p' | wc -l | tr -d ' ')"
chk "it is not flagged SOMATIC" 0 "$(echo "$rec2" | cut -f8 | grep -c 'SOMATIC')"

# ------------------------------------------------- shared phasing (M4) -----
read -r n_ps n_flip <<<"$(bcftools query -f '[%GT\t%PS\t]\n' "$T/paired.vcf" 2>/dev/null | awk -F'\t' '
  $2!="." && $4!="." { if ($2!=$4) m++
    if (($1=="0|1"||$1=="1|0") && ($3=="0|1"||$3=="1|0")) { if ($1!=$3) f++ } }
  END{print m+0, f+0}')"
chk "no phase-set disagreement between samples" 0 "$n_ps"
chk "no orientation flip at shared heterozygous sites" 0 "$n_flip"

echo
if [[ $fail -eq 0 ]]; then echo "  paired VCF: all $pass assertions passed"; exit 0
else echo "  paired VCF: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
