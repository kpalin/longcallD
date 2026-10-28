#!/usr/bin/env bash
#
# QUAL, GQ and PL (src/var_qual.c).
#
# The quality model only scores records, so most of what can be checked on the
# bundled data is internal consistency -- PL against GQ against the printed
# genotype, UPR against the slot grid -- and direction: more reads, a cleaner
# repeat context and both strands must each score higher than their opposite.
# Whether the numbers are calibrated is a benchmark question; see
# eval/qual_eval.py and the README.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
REF=$ROOT/test_data/chr11_2M.fa
HIFI=$ROOT/test_data/HG002_chr11_hifi_test.bam
T=${KEEP_T:-$(mktemp -d)}; [[ -z ${KEEP_T:-} ]] && trap 'rm -rf "$T"' EXIT; mkdir -p "$T"
HET=1258256; HET_ALT=G; SITE=1258300; NEW=A; PVAL=0.01
HPHET=1269349; HPHET_ALT=T; HPSITE=1268705; CLEANSITE=1269500

pass=0; fail=0
chk() {
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

# Every sample column of every record: GQ must be what its own PL and its own
# printed genotype give -- the best other genotype's PL minus the printed one's,
# floored at 0 and capped at 99 -- and the smallest PL must be 0.
gq_check() { # $1 = vcf; prints the number of sample columns that disagree
  awk -F'\t' '!/^#/ {
    n = split($9, k, ":"); for (i = 1; i <= n; i++) ix[k[i]] = i
    for (s = 10; s <= NF; s++) {
      split($s, v, ":"); gt = v[ix["GT"]]; if (gt ~ /\./) continue
      d = (substr(gt,1,1) != "0") + (substr(gt,3,1) != "0")
      split(v[ix["PL"]], pl, ","); mn = pl[1]; if (pl[2] < mn) mn = pl[2]; if (pl[3] < mn) mn = pl[3]
      o = 1e9; for (g = 0; g <= 2; g++) if (g != d && pl[g+1] < o) o = pl[g+1]
      want = o - pl[d+1]; if (want < 0) want = 0; if (want > 99) want = 99
      if (mn != 0 || v[ix["GQ"]] != want) bad++
    }
  } END { print bad + 0 }' "$1"
}

# ---------------------------------------------------------- single sample ---
# --qual-weights none: the likelihood alone, which is what the depth and QLIK
# checks below are about. The fitted model is exercised further down.
"$BIN" call --hifi --qual-weights none "$REF" "$HIFI" -t 1 -o "$T/full.vcf" 2>/dev/null
samtools view -b -s 7.5 "$HIFI" -o "$T/half.bam" && samtools index "$T/half.bam"
"$BIN" call --hifi --qual-weights none "$REF" "$T/half.bam" -t 1 -o "$T/half.vcf" 2>/dev/null
nrec=$(grep -vc '^#' "$T/full.vcf")
chk "single: records are written" yes "$([[ $nrec -gt 100 ]] && echo yes || echo no)"
for tag in QLIK SOR MQRS BQRS RU RPL; do
  chk "INFO $tag is declared" 1 "$(grep -c "##INFO=<ID=$tag," "$T/full.vcf")"
done
for tag in PL UPR UPF; do
  chk "FORMAT $tag is declared" 1 "$(grep -c "##FORMAT=<ID=$tag," "$T/full.vcf")"
done
chk "bcftools accepts the records" 0 "$(bcftools view "$T/full.vcf" >/dev/null 2>&1; echo $?)"
chk "single: GQ follows from PL and GT in every record" 0 "$(gq_check "$T/full.vcf")"
# With no model set, QUAL is QLIK to one decimal (and capped at 999, where QLIK is not).
chk "single: without a model QUAL is QLIK" 0 \
    "$(awk -F'\t' '!/^#/ { if (match($8, /QLIK=[0-9.]+/)) { q = substr($8, RSTART+5, RLENGTH-5) + 0; if (q > 999) q = 999; if ($6 - q > 0.06 || q - $6 > 0.06) b++ } } END {print b+0}' "$T/full.vcf")"
# The built-in model must still produce valid records.
"$BIN" call --hifi "$REF" "$HIFI" -t 1 -o "$T/model.vcf" 2>/dev/null
chk "single: the built-in model scores the same records" "$nrec" "$(grep -vc '^#' "$T/model.vcf")"
chk "single: the built-in model keeps GQ consistent with PL" 0 "$(gq_check "$T/model.vcf")"
# Halving the reads must not raise QUAL on the sites both runs call, and must
# lower it on most of them: sampling error is what QLIK measures.
cmp_depth=$(awk -F'\t' 'NR==FNR { if (!/^#/) q[$2" "$4" "$5] = $6; next }
  !/^#/ && ($2" "$4" "$5) in q { n++; if (q[$2" "$4" "$5] > $6) lo++; if (q[$2" "$4" "$5] < $6) hi++ }
  END { printf "%d %d %d", n, lo+0, hi+0 }' "$T/full.vcf" "$T/half.vcf")
read -r n_shared n_lower n_higher <<< "$cmp_depth"
chk "depth: half the reads scores lower on most shared sites" yes \
    "$([[ $n_shared -gt 50 && $((n_lower * 10)) -ge $((n_shared * 8)) ]] && echo yes || echo "no ($cmp_depth)")"
chk "depth: and higher on almost none" yes \
    "$([[ $((n_higher * 20)) -le $n_shared ]] && echo yes || echo "no ($cmp_depth)")"
# The unphased share must count reads, not invent them.
chk "single: UPR never exceeds DP" 0 \
    "$(awk -F'\t' '!/^#/ { n = split($9, k, ":"); for (i = 1; i <= n; i++) ix[k[i]] = i
        split($10, v, ":"); split(v[ix["UPR"]], u, ","); if (u[1] > v[ix["DP"]] || u[2] > u[1]) b++ } END {print b+0}' "$T/full.vcf")"
chk "--min-qual 1000 marks every record LowQual" "$nrec" \
    "$("$BIN" call --hifi --min-qual 1000 "$REF" "$HIFI" -t 1 2>/dev/null | awk -F'\t' '!/^#/ && $7=="LowQual"' | wc -l | tr -d ' ')"
chk "--min-qual off by default" 0 "$(awk -F'\t' '!/^#/ && $7=="LowQual"' "$T/full.vcf" | wc -l | tr -d ' ')"

# ----------------------------------------------------------------- paired ---
samtools view "$HIFI" | awk '{print $1}' | sort -u > "$T/names"
n=$(wc -l < "$T/names"); head -n $((n/2)) "$T/names" > "$T/n1"; tail -n +$((n/2+1)) "$T/names" > "$T/n2"
subset() { { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$1"; }; }
subset "$T/n2" | samtools sort -o "$T/normal.bam" - 2>/dev/null; samtools index "$T/normal.bam"
subset "$T/n1" | python3 -I "$ROOT/eval/inject_somatic.py" "$HET" "$HET_ALT" "$SITE" "$NEW" 1.0 2>/dev/null \
  > "$T/tumor.sam"
samtools sort -o "$T/tumor.bam" "$T/tumor.sam" 2>/dev/null; samtools index "$T/tumor.bam"
# The same tumour with every reverse-strand read that carries the injected base
# removed: the alt is then seen on one strand only, the reference on both.
cat > "$T/fwd.py" <<'PY'
import sys
sys.dont_write_bytecode = True  # no __pycache__ in eval/
sys.path.insert(0, sys.argv[1])
from inject_somatic import CIG, query_offset
site, new = int(sys.argv[2]), sys.argv[3]
for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line); continue
    f = line.split("\t")
    if int(f[1]) & 16:
        q = query_offset(int(f[3]), [(int(n), op) for n, op in CIG.findall(f[5])], site)
        if q is not None and f[9][q] == new:
            continue
    sys.stdout.write(line)
PY
python3 -I "$T/fwd.py" "$ROOT/eval" "$SITE" "$NEW" < "$T/tumor.sam" > "$T/tumor_fwd.sam"
samtools sort -o "$T/tumor_fwd.bam" "$T/tumor_fwd.sam" 2>/dev/null; samtools index "$T/tumor_fwd.bam"

"$BIN" call --hifi --som-pval $PVAL "$REF" "$T/tumor.bam" -N "$T/normal.bam" -t 1 -o "$T/paired.vcf" 2>/dev/null
"$BIN" call --hifi --som-pval $PVAL "$REF" "$T/tumor_fwd.bam" -N "$T/normal.bam" -t 1 -o "$T/paired_fwd.vcf" 2>/dev/null
info_at() { awk -F'\t' -v p="$2" -v k="$3" '!/^#/ && $2==p { n = split($8, a, ";"); for (i = 1; i <= n; i++) if (index(a[i], k"=") == 1) print substr(a[i], length(k)+2) }' "$1"; }
qual_at() { awk -F'\t' -v p="$2" '!/^#/ && $2==p {print $6}' "$1"; }

chk "paired: GQ follows from each sample's own PL and GT" 0 "$(gq_check "$T/paired.vcf")"
chk "paired: the injected SNV is somatic" yes \
    "$(awk -F'\t' -v p=$SITE '!/^#/ && $2==p && $8 ~ /(^|;)SOMATIC(;|$)/' "$T/paired.vcf" | grep -q . && echo yes || echo no)"
chk "paired: a clonal somatic SNV against a clean 14x normal has QLIK >= 20" yes \
    "$(awk -v q="$(info_at "$T/paired.vcf" $SITE QLIK)" 'BEGIN{print (q >= 20) ? "yes" : "no (" q ")"}')"
# The normal is 0/0 and the tumour 0/1 there; each GQ must speak for its own column.
chk "paired: the two samples get their own GQ" yes \
    "$(awk -F'\t' -v p=$SITE '!/^#/ && $2==p { n = split($9, k, ":"); for (i = 1; i <= n; i++) ix[k[i]] = i
        split($10, a, ":"); split($11, b, ":"); print (a[ix["PL"]] != b[ix["PL"]]) ? "yes" : "no" }' "$T/paired.vcf")"
# UPR[0] counts the reads the grid left in its unphased slot, so it is DP - HDP.
chk "paired: UPR equals DP minus the two haplotype depths" 0 \
    "$(awk -F'\t' '!/^#/ { n = split($9, k, ":"); for (i = 1; i <= n; i++) ix[k[i]] = i
        for (s = 10; s <= 11; s++) { split($s, v, ":"); split(v[ix["HDP"]], h, ","); split(v[ix["UPR"]], u, ",")
          if (u[1] != v[ix["DP"]] - h[1] - h[2]) b++ } } END {print b+0}' "$T/paired.vcf")"
sor_both=$(info_at "$T/paired.vcf" $SITE SOR); sor_fwd=$(info_at "$T/paired_fwd.vcf" $SITE SOR)
chk "strand: SOR rises when the alt is on one strand" yes \
    "$(awk -v a="$sor_both" -v b="$sor_fwd" 'BEGIN{print (b != "" && b > a + 1) ? "yes" : "no (" a " -> " b ")"}')"
# A strand weight must act on QUAL through --qual-weights, and only there. Both
# classes get it: dropping the reverse-strand alt reads leaves too few for the
# classifier, which can then label the site Ambiguous and score it as germline.
printf 'hifi %s bias -6\nhifi %s one_strand 5\nhifi %s sor 1\n' som_snv som_snv som_snv germ_snv germ_snv germ_snv > "$T/w.txt"
q_both=$("$BIN" call --hifi --som-pval $PVAL --qual-weights "$T/w.txt" "$REF" "$T/tumor.bam" -N "$T/normal.bam" -t 1 2>/dev/null \
         | awk -F'\t' -v p=$SITE '!/^#/ && $2==p {print $6}')
q_fwd=$("$BIN" call --hifi --som-pval $PVAL --qual-weights "$T/w.txt" "$REF" "$T/tumor_fwd.bam" -N "$T/normal.bam" -t 1 2>/dev/null \
        | awk -F'\t' -v p=$SITE '!/^#/ && $2==p {print $6}')
chk "strand: with a strand weight, one-strand support scores lower" yes \
    "$(awk -v a="$q_both" -v b="$q_fwd" 'BEGIN{print (a != "" && b != "" && b + 0 < a + 0) ? "yes" : "no (" a " -> " b ")"}')"
chk "a malformed --qual-weights file is refused" yes \
    "$(printf 'hifi som_snv nosuchfeature 1\n' > "$T/bad.txt"; "$BIN" call --hifi --qual-weights "$T/bad.txt" "$REF" "$HIFI" -t 1 >/dev/null 2>&1 && echo no || echo yes)"

# ---------------------------------------------------------- repeat context ---
del_build() { # $1 = site, $2 = output bam
  subset "$T/n1" | python3 -I "$ROOT/eval/inject_somatic.py" "$HPHET" "$HPHET_ALT" "$1" "del:1" 1.0 2>/dev/null \
    | samtools sort -o "$2" - 2>/dev/null
  samtools index "$2"
}
del_build "$HPSITE" "$T/tumor_hp.bam"; del_build "$CLEANSITE" "$T/tumor_ci.bam"
# Artifact mask off, so both deletions are SOMATIC and differ only in context.
"$BIN" call --hifi --som-pval $PVAL --som-artifact-mask-indel 0 "$REF" "$T/tumor_hp.bam" -N "$T/normal.bam" -t 1 -o "$T/hp.vcf" 2>/dev/null
"$BIN" call --hifi --som-pval $PVAL --som-artifact-mask-indel 0 "$REF" "$T/tumor_ci.bam" -N "$T/normal.bam" -t 1 -o "$T/ci.vcf" 2>/dev/null
hp_anchor=$((HPSITE-1)); ci_anchor=$((CLEANSITE-1))
chk "repeat: the homopolymer deletion reports its tract" "T 9" \
    "$(info_at "$T/hp.vcf" $hp_anchor RU) $(info_at "$T/hp.vcf" $hp_anchor RPL)"
chk "repeat: the clean-sequence deletion reports none" "" "$(info_at "$T/ci.vcf" $ci_anchor RU)"
ql_hp=$(info_at "$T/hp.vcf" $hp_anchor QLIK); ql_ci=$(info_at "$T/ci.vcf" $ci_anchor QLIK)
chk "repeat: a deletion in a 9-base homopolymer scores lower than one in clean sequence" yes \
    "$(awk -v a="$ql_hp" -v b="$ql_ci" 'BEGIN{print (a != "" && b != "" && a + 0 < b + 0) ? "yes" : "no (" a " vs " b ")"}')"

echo
if [[ $fail -eq 0 ]]; then echo "  quality: all $pass assertions passed"; exit 0
else echo "  quality: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
