#!/usr/bin/env bash
#
# Invariants the phased alignment output must satisfy. Run after any change to
# the BAM writer.
#
# Baseline before the reorder buffer existed: --refine-aln on the bundled ONT
# file emitted 6 of 587 records out of coordinate order while @HD still said
# SO:coordinate, so samtools index failed with
#   [E::hts_idx_push] Unsorted positions on sequence #11: 1308835 followed by 1308804
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$ROOT/bin/longcallD}
REF=$ROOT/test_data/chr11_2M.fa
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

pass=0; fail=0
chk() {
  if [[ "$2" == "$3" ]]; then printf '  PASS  %s\n' "$1"; pass=$((pass+1))
  else printf '  FAIL  %s (got %q, want %q)\n' "$1" "$3" "$2"; fail=$((fail+1)); fi
}
# records that start before their predecessor on the same reference
n_unsorted() { samtools view "$1" | awk '{if($3==p3 && $4<p4) n++; p3=$3; p4=$4} END{print n+0}'; }

[[ -x $BIN ]] || { echo "no binary at $BIN; run make first" >&2; exit 1; }

for rt in hifi ont; do
  flag=""; [[ $rt == ont ]] && flag="--ont"
  BAM=$ROOT/test_data/HG002_chr11_${rt}_test.bam
  n_in=$(samtools view -c "$BAM")

  for t in 1 4; do
    for refine in "" "--refine-aln"; do
      tag="$rt -t$t${refine:+ --refine-aln}"
      out="$T/${rt}_t${t}${refine:+_r}.bam"
      "$BIN" call $flag $refine "$REF" "$BAM" -t "$t" -b "$out" -o /dev/null >/dev/null 2>&1
      chk "$tag writes a BAM" yes "$([[ -s $out ]] && echo yes || echo no)"
      chk "$tag coordinate-sorted" 0 "$(n_unsorted "$out")"
      chk "$tag header says SO:coordinate" "SO:coordinate" \
          "$(samtools view --no-PG -H "$out" | grep -oE 'SO:[a-z]+' | head -1)"
      rm -f "$out.bai"
      samtools index "$out" 2>/dev/null
      chk "$tag samtools index succeeds" yes "$([[ -f $out.bai ]] && echo yes || echo no)"
      chk "$tag conserves record count" "$n_in" "$(samtools view -c "$out")"
    done
  done

  # --no-reorder must not claim to be sorted.
  out=$T/${rt}_nr.bam
  "$BIN" call $flag --refine-aln --no-reorder "$REF" "$BAM" -t 1 -b "$out" -o /dev/null >/dev/null 2>&1
  chk "$rt --no-reorder header says SO:unsorted" "SO:unsorted" \
      "$(samtools view --no-PG -H "$out" | grep -oE 'SO:[a-z]+' | head -1)"
  chk "$rt --no-reorder conserves record count" "$n_in" "$(samtools view -c "$out")"

  # Reordering must not alter any record, only their order.
  chk "$rt reordered output has the same records as --no-reorder" same \
      "$([[ "$(samtools view "$T/${rt}_t1_r.bam" | sort | md5sum)" == "$(samtools view "$out" | sort | md5sum)" ]] && echo same || echo differs)"

  # -t 1 and -t 4 must agree.
  chk "$rt deterministic across thread counts" same \
      "$([[ "$(samtools view "$T/${rt}_t1_r.bam" | md5sum)" == "$(samtools view "$T/${rt}_t4_r.bam" | md5sum)" ]] && echo same || echo differs)"
done

echo
if [[ $fail -eq 0 ]]; then echo "  BAM output: all $pass assertions passed"; exit 0
else echo "  BAM output: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
