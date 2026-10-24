#!/usr/bin/env bash
#
# Behavioural assertions for the option parser. Every case here is a bug that
# was reproducible against the built binary before the fix, so this doubles as
# a regression test and as documentation of what each option is supposed to do.
#
# Runs in a scratch directory because one of the bugs was an option *creating a
# file named after its own argument*, which is only detectable by looking at the
# filesystem afterwards.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
REF=$ROOT/test_data/chr11_2M.fa
BAM=$ROOT/test_data/HG002_chr11_hifi_test.bam
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"

pass=0; fail=0
chk() { # name expected actual
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
nvar() { grep -vc '^#' "$1" 2>/dev/null || echo 0; }

[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

# --max-cov shared 'C' with --out-cram: it used to open a CRAM named "1000".
"$BIN" call --max-cov 1000 "$REF" "$BAM" -o mc.vcf >/dev/null 2>&1
chk "--max-cov exits 0"                       0         "$?"
chk "--max-cov creates no file named '1000'"  absent    "$([[ -e 1000 ]] && echo present || echo absent)"
chk "--max-cov still writes a VCF"            yes       "$([[ -s mc.vcf ]] && echo yes || echo no)"

# ...and it must actually take effect, not merely be accepted.
"$BIN" call --max-cov 5 "$REF" "$BAM" -o mc5.vcf >/dev/null 2>&1
chk "--max-cov changes the call set" differs \
    "$([[ "$(nvar mc.vcf)" != "$(nvar mc5.vcf)" ]] && echo differs || echo same)"

# These three exited 0 with no VCF, or silently set a different option.
for o in "--clip-flank 100" "--merge-dis 100" "--end-clip 150"; do
  "$BIN" call $o "$REF" "$BAM" -o t.vcf >/dev/null 2>&1; ec=$?
  chk "$o works" "0:yes" "$ec:$([[ $(nvar t.vcf) -gt 0 ]] && echo yes || echo no)"
done

# Retired options: backing fields no longer exist, so they must be rejected.
for o in "--noisy-flank 10" "--full-read 5"; do
  "$BIN" call $o "$REF" "$BAM" -o r.vcf >/dev/null 2>&1
  chk "$o rejected non-zero" 1 "$?"
done

# --noisy-rat was retired with the rest, then brought back with a real field
# when the re-alignment trigger was separated from -a. It must work, and it
# must not be a synonym for -a: a region-noise threshold of 1.0 re-aligns
# nothing, while -a 1.0 would also call nothing.
"$BIN" call --noisy-rat 0.3 "$REF" "$BAM" -o t.vcf >/dev/null 2>&1; ec=$?
chk "--noisy-rat 0.3 works" "0:yes" "$ec:$([[ $(nvar t.vcf) -gt 0 ]] && echo yes || echo no)"
"$BIN" call --noisy-rat 1.0 "$REF" "$BAM" -o t.vcf >/dev/null 2>&1
chk "--noisy-rat 1.0 still calls variants" yes "$([[ $(nvar t.vcf) -gt 0 ]] && echo yes || echo no)"

# A typo must not look like a successful run that happened to produce nothing.
"$BIN" call --no-such-option "$REF" "$BAM" >/dev/null 2>&1
chk "unknown option exits non-zero" 1 "$?"

# case 'g' fell through to case 't', setting n_threads = atoi("left") = 0.
"$BIN" call --gap-aln left "$REF" "$BAM" -o g.vcf >/dev/null 2>&1; ec=$?
chk "--gap-aln left works" "0:yes" "$ec:$([[ -s g.vcf ]] && echo yes || echo no)"
"$BIN" call --gap-aln bogus "$REF" "$BAM" -o g2.vcf >/dev/null 2>&1
chk "--gap-aln bogus rejected" 1 "$?"

# opt_str said "L:" while the table said no_argument, so -L ate ref.fa.
printf '%s\n' "$BAM" > list.txt
"$BIN" call -L "$REF" list.txt -o l.vcf >/dev/null 2>&1; ec=$?
chk "-L ref.fa list.txt works" "0:yes" "$ec:$([[ -s l.vcf ]] && echo yes || echo no)"

# sam_hdr_add_pg() added one @PG per PP-chain leaf: 87 on the ONT test file.
for rt in hifi ont; do
  flag=""; [[ $rt == ont ]] && flag="--ont"
  "$BIN" call $flag "$REF" "$ROOT/test_data/HG002_chr11_${rt}_test.bam" \
      -b "pg_$rt.bam" -o /dev/null >/dev/null 2>&1
  chk "$rt: exactly one PN:longcallD @PG" 1 \
      "$(samtools view --no-PG -H "pg_$rt.bam" | grep -c 'PN:longcallD')"
  chk "$rt: output has input @PG count + 1" \
      "$(( $(samtools view --no-PG -H "$ROOT/test_data/HG002_chr11_${rt}_test.bam" | grep -c '^@PG') + 1 ))" \
      "$(samtools view --no-PG -H "pg_$rt.bam" | grep -c '^@PG')"
done

echo
if [[ $fail -eq 0 ]]; then echo "  CLI: all $pass assertions passed"; exit 0
else echo "  CLI: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
