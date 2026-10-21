#!/usr/bin/env bash
#
# Per-read sample origin: every read the caller uses must be attributed to the
# input file it came from.
#
# The -N files are appended to in_bam_fns only after set_bam_list_file() has
# finished permuting that array, so sample identity reduces to
# "index >= n_tumor_bam_fn". These assertions exist because that invariant is
# easy to break: an off-by-one in the boundary silently mixes tumour and normal
# evidence, which would not show up as a crash, only as wrong calls.
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

# Split the bundled file into disjoint read-name sets, so the two "samples"
# share no molecule and the expected per-sample tallies are unambiguous.
samtools view "$HIFI" | awk '{print $1}' | sort -u > "$T/names"
n=$(wc -l < "$T/names")
split_at() { # $1=count for the first part, $2/$3 = output bam names
  head -n "$1" "$T/names" > "$T/n1"; tail -n +"$(($1+1))" "$T/names" > "$T/n2"
  local i=0
  for part in n1 n2; do
    i=$((i+1)); local out="$T/$([[ $i == 1 ]] && echo "$2" || echo "$3")"
    { samtools view --no-PG -H "$HIFI"; samtools view "$HIFI" | grep -Fwf "$T/$part"; } \
      | samtools sort -o "$out" - 2>/dev/null
    samtools index "$out"
  done
}
# Reads the caller actually uses: primary only (-F 2308 drops unmapped,
# secondary and supplementary) at or above the default MAPQ floor.
callable() { samtools view -c -F 2308 -q 30 "$1"; }
tally() { # $1=tumour bam $2..=normal bams -> "tumour N, normal M"
  local t=$1; shift
  local args=(); for f in "$@"; do args+=(-N "$f"); done
  "$BIN" call "$REF" "$t" "${args[@]}" -t 1 -o /dev/null 2>&1 \
    | grep -oE 'Reads used for calling: .*' | sed 's/Reads used for calling: //'
}

split_at $((n/2)) half_a.bam half_b.bam
chk "even split attributes reads correctly" \
    "tumour $(callable "$T/half_a.bam"), normal $(callable "$T/half_b.bam")" \
    "$(tally "$T/half_a.bam" "$T/half_b.bam")"

split_at $((n/4)) quarter.bam rest.bam
chk "asymmetric split attributes reads correctly" \
    "tumour $(callable "$T/quarter.bam"), normal $(callable "$T/rest.bam")" \
    "$(tally "$T/quarter.bam" "$T/rest.bam")"
chk "roles swap cleanly" \
    "tumour $(callable "$T/rest.bam"), normal $(callable "$T/quarter.bam")" \
    "$(tally "$T/rest.bam" "$T/quarter.bam")"

split_at $((n/2)) half_a.bam half_b.bam
chk "-N twice sums both normals" \
    "tumour $(callable "$T/quarter.bam"), normal $(( $(callable "$T/half_a.bam") + $(callable "$T/half_b.bam") ))" \
    "$(tally "$T/quarter.bam" "$T/half_a.bam" "$T/half_b.bam")"

# Without -N nothing paired should be announced at all.
chk "single sample stays single-sample" 0 \
    "$("$BIN" call "$REF" "$HIFI" -t 1 -o /dev/null 2>&1 | grep -cE 'Paired input|Reads used for calling')"

# ---------------------------------------------------------------------------
# Multi-input alignment output must stay coordinate-sorted.
#
# The writer used to iterate one input file to exhaustion before starting the
# next, which groups the output by file: with these two halves as -X inputs the
# position jumped backwards 185 kb at the boundary, while @HD still claimed
# SO:coordinate. Records are now emitted by a k-way merge across the files.
# ---------------------------------------------------------------------------
n_unsorted() { samtools view "$1" | awk '{if($3==p3 && $4<p4) n++; p3=$3; p4=$4} END{print n+0}'; }
n_half_phased() { samtools view "$1" | awk '{hp=/HP:i:/?1:0; ps=/PS:i:/?1:0; if(hp!=ps) n++} END{print n+0}'; }
n_both=$(( $(samtools view -c "$T/half_a.bam") + $(samtools view -c "$T/half_b.bam") ))

# -X pools both files as one sample, so -b holds everything. -N makes the
# second file a separate sample, so -b holds only the tumour and the normal
# goes to --out-normal-bam (covered further down).
n_tumor_only=$(samtools view -c "$T/half_a.bam")
for mode in -X -N; do
  for refine in "" "--refine-aln"; do
    tag="$mode${refine:+ --refine-aln}"
    out=$T/multi.bam
    expect=$([[ $mode == -X ]] && echo "$n_both" || echo "$n_tumor_only")
    "$BIN" call $refine "$REF" "$T/half_a.bam" $mode "$T/half_b.bam" -t 1 -b "$out" -o /dev/null >/dev/null 2>&1
    chk "$tag -b holds the expected records" "$expect" "$(samtools view -c "$out")"
    chk "$tag output is coordinate-sorted" 0 "$(n_unsorted "$out")"
    chk "$tag output keeps HP/PS joint" 0 "$(n_half_phased "$out")"
    rm -f "${out:?}.bai"; samtools index "$out" 2>/dev/null
    chk "$tag output indexes" yes "$([[ -f $out.bai ]] && echo yes || echo no)"
  done
done

# ---------------------------------------------------------------------------
# Phasing is anchored only on variants that look germline heterozygous in the
# NORMAL sample.
#
# That is what makes haplotype index h denote the same parental chromosome in
# both samples: a tumour-only variant would otherwise define a "haplotype" that
# exists only in the tumour, and comparing slot (T,h) against slot (N,h) would
# be comparing two different things. Variants are still *called* in both
# samples; this only governs what the phasing is anchored on.
#
# Driving --min-normal-phase-dp to an unreachable value makes every variant
# ineligible, so phasing must collapse while calling continues. If phasing
# survived that, the filter would not be wired into it at all.
# ---------------------------------------------------------------------------
n_phased() { samtools view "$1" | awk '{if(/HP:i:/ && /PS:i:/) n++} END{print n+0}'; }

"$BIN" call "$REF" "$T/half_a.bam" -N "$T/half_b.bam" -t 1 -b "$T/ph_def.bam" -o "$T/ph_def.vcf" >/dev/null 2>&1
chk "paired: phasing happens with the default threshold" yes     "$([[ $(n_phased "$T/ph_def.bam") -gt 0 ]] && echo yes || echo no)"
chk "paired: variants are called with the default threshold" yes     "$([[ $(grep -vc '^#' "$T/ph_def.vcf") -gt 0 ]] && echo yes || echo no)"

"$BIN" call "$REF" "$T/half_a.bam" -N "$T/half_b.bam" -t 1 --min-normal-phase-dp 100000 \
    -b "$T/ph_none.bam" -o "$T/ph_none.vcf" >/dev/null 2>&1
chk "paired: no eligible variant means no phasing" 0 "$(n_phased "$T/ph_none.bam")"
chk "paired: calling survives with nothing to phase on" yes     "$([[ $(grep -vc '^#' "$T/ph_none.vcf") -gt 0 ]] && echo yes || echo no)"
chk "paired: no eligible variant still exits cleanly" 0     "$("$BIN" call "$REF" "$T/half_a.bam" -N "$T/half_b.bam" -t 1 --min-normal-phase-dp 100000 -o /dev/null >/dev/null 2>&1; echo $?)"

# The option must not touch single-sample behaviour at all.
"$BIN" call "$REF" "$HIFI" -t 1 --min-normal-phase-dp 100000 -b "$T/ss_opt.bam" -o /dev/null >/dev/null 2>&1
"$BIN" call "$REF" "$HIFI" -t 1 -b "$T/ss_noopt.bam" -o /dev/null >/dev/null 2>&1
chk "single sample ignores --min-normal-phase-dp" same \
    "$([[ "$(samtools view "$T/ss_opt.bam" | md5sum)" == "$(samtools view "$T/ss_noopt.bam" | md5sum)" ]] && echo same || echo differs)"

# ---------------------------------------------------------------------------
# Separate output BAMs for tumour and normal.
#
# The requirement is that phasing is consistent between them, so the same phase
# sets must appear in both files. The risk is the opposite: global_read_i
# indexes chunk->haps across every input, so routing records to two files is
# exactly where an off-by-N puts one read's haplotype on another read. The
# per-read signature check would abort if that happened.
# ---------------------------------------------------------------------------
n_t=$(samtools view -c "$T/half_a.bam"); n_n=$(samtools view -c "$T/half_b.bam")
"$BIN" call --refine-aln "$REF" "$T/half_a.bam" -N "$T/half_b.bam" -t 1 \
    -b "$T/o_t.bam" --out-normal-bam "$T/o_n.bam" -o /dev/null >/dev/null 2>&1
chk "two output BAMs: tumour gets exactly its own records" "$n_t" "$(samtools view -c "$T/o_t.bam")"
chk "two output BAMs: normal gets exactly its own records" "$n_n" "$(samtools view -c "$T/o_n.bam")"
chk "two output BAMs: no read appears in both" 0 \
    "$(comm -12 <(samtools view "$T/o_t.bam" | cut -f1 | sort -u) <(samtools view "$T/o_n.bam" | cut -f1 | sort -u) | wc -l | tr -d ' ')"
for f in o_t o_n; do
  chk "$f is coordinate-sorted" 0 "$(n_unsorted "$T/$f.bam")"
  chk "$f keeps HP/PS joint" 0 "$(n_half_phased "$T/$f.bam")"
  rm -f "${T:?}/${f:?}.bam.bai"; samtools index "$T/$f.bam" 2>/dev/null
  chk "$f indexes" yes "$([[ -f $T/$f.bam.bai ]] && echo yes || echo no)"
done
# The stated requirement: phasing consistent between the two output BAMs.
chk "both output BAMs use the same phase sets" same \
    "$([[ "$(samtools view "$T/o_t.bam" | grep -oE 'PS:i:[0-9]+' | sort -u | md5sum)" == \
          "$(samtools view "$T/o_n.bam" | grep -oE 'PS:i:[0-9]+' | sort -u | md5sum)" ]] && echo same || echo differs)"

# With only -b, the normal's reads must not leak into the tumour's file -- and
# read_i must still advance over them, which the signature check enforces.
"$BIN" call --refine-aln "$REF" "$T/half_a.bam" -N "$T/half_b.bam" -t 1 -b "$T/only_t.bam" -o /dev/null >/dev/null 2>&1
chk "with only -b, output holds just the tumour" "$n_t" "$(samtools view -c "$T/only_t.bam")"
chk "with only -b, no normal read leaks in" 0 \
    "$(comm -12 <(samtools view "$T/only_t.bam" | cut -f1 | sort -u) <(samtools view "$T/half_b.bam" | cut -f1 | sort -u) | wc -l | tr -d ' ')"

chk "--out-normal-bam without -N is refused" 1 \
    "$("$BIN" call "$REF" "$HIFI" -t 1 --out-normal-bam "$T/nope.bam" -o /dev/null >/dev/null 2>&1; echo $?)"

echo
if [[ $fail -eq 0 ]]; then echo "  paired input: all $pass assertions passed"; exit 0
else echo "  paired input: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
