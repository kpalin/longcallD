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
# records carrying exactly one of HP/PS -- the pair is meaningless alone
n_half_phased() { samtools view "$1" | awk '{hp=/HP:i:/?1:0; ps=/PS:i:/?1:0; if(hp!=ps) n++} END{print n+0}'; }

# tagged supplementary records whose HP/PS disagree with their own primary
n_supp_disagree() {
  samtools view "$1" | awk '
    { hp=""; ps=""
      for(i=12;i<=NF;i++){ if($i ~ /^HP:i:/){split($i,a,":"); hp=a[3]}
                           if($i ~ /^PS:i:/){split($i,b,":"); ps=b[3]} }
      if (and($2,2048)) { s_hp[$1]=hp; s_ps[$1]=ps }
      else if (!and($2,256) && !and($2,4)) { p_hp[$1]=hp; p_ps[$1]=ps }
    }
    END { for (k in s_hp) if (s_hp[k]!="" && (s_hp[k]!=p_hp[k] || s_ps[k]!=p_ps[k])) n++
          print n+0 }'
}

# secondary or unmapped records carrying our tags (they never should)
n_sec_unmapped_tagged() {
  samtools view "$1" | awk '(and($2,256) || and($2,4)) && (/HP:i:/ || /PS:i:/) {n++} END{print n+0}'
}

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
      # Baseline before the fix: 1 HP-only + 5 PS-only on ONT, both read types
      # affected equally with and without --refine-aln.
      chk "$tag every record has both HP and PS or neither" 0 "$(n_half_phased "$out")"
      # Supplementary alignments inherit their read's primary phase; a
      # supplementary whose primary is in another chunk gets nothing.
      chk "$tag supplementary tags agree with their primary" 0 "$(n_supp_disagree "$out")"
      chk "$tag secondary/unmapped records carry no HP/PS" 0 "$(n_sec_unmapped_tagged "$out")"
    done
  done

  # --no-reorder must not claim to be sorted.
  out=$T/${rt}_nr.bam
  "$BIN" call $flag --refine-aln --no-reorder "$REF" "$BAM" -t 1 -b "$out" -o /dev/null >/dev/null 2>&1
  chk "$rt --no-reorder header says SO:unsorted" "SO:unsorted" \
      "$(samtools view --no-PG -H "$out" | grep -oE 'SO:[a-z]+' | head -1)"
  chk "$rt --no-reorder conserves record count" "$n_in" "$(samtools view -c "$out")"
  chk "$rt --no-reorder keeps HP/PS joint" 0 "$(n_half_phased "$out")"

  # Reordering must not alter any record, only their order.
  chk "$rt reordered output has the same records as --no-reorder" same \
      "$([[ "$(samtools view "$T/${rt}_t1_r.bam" | sort | md5sum)" == "$(samtools view "$out" | sort | md5sum)" ]] && echo same || echo differs)"

  # -t 1 and -t 4 must agree.
  chk "$rt deterministic across thread counts" same \
      "$([[ "$(samtools view "$T/${rt}_t1_r.bam" | md5sum)" == "$(samtools view "$T/${rt}_t4_r.bam" | md5sum)" ]] && echo same || echo differs)"
done

# ---------------------------------------------------------------------------
# Aux-tag integrity after CIGAR refinement.
#
# refine_bam1() performs hand-rolled surgery on b->data, and update_bam1_tags()
# rewrites NM/MD/cs from the new alignment. Three things must hold: records whose
# alignment did not change must come out byte-identical, only the derived tags may
# change on records that did, and the rewritten values must be right -- checked
# against samtools calmd, an independent implementation.
# ---------------------------------------------------------------------------
for rt in hifi ont; do
  flag=""; [[ $rt == ont ]] && flag="--ont"
  BAM=$ROOT/test_data/HG002_chr11_${rt}_test.bam
  out=$T/${rt}_t1_r.bam   # already produced above, with --refine-aln at -t 1

  # Pair by qname+flag and compare POS, CIGAR and the aux fields other than ours.
  aux_of() { samtools view "$1" | awk '{a=""; for(i=12;i<=NF;i++) if($i !~ /^(HP|PS):i:/) a=a" "$i; print $1"_"$2"\t"$4"\t"$6"\t"a}' | sort; }
  join -t$'\t' -j1 <(aux_of "$BAM") <(aux_of "$out") > "$T/${rt}.pair"

  chk "$rt untouched alignments keep every aux tag byte-identical" 0 \
      "$(awk -F'\t' '$2==$5 && $3==$6 && $4!=$7 {n++} END{print n+0}' "$T/${rt}.pair")"

  # Only NM, MD and cs may differ where the alignment changed.
  chk "$rt refined records change only NM/MD/cs" "" \
      "$(awk -F'\t' '$2!=$5 || $3!=$6 {
            n=split($4,a," "); m=split($7,b," "); delete iv; delete ov
            for(i=1;i<=n;i++){split(a[i],x,":"); iv[x[1]]=a[i]}
            for(i=1;i<=m;i++){split(b[i],y,":"); ov[y[1]]=b[i]}
            for(k in iv) if (iv[k]!=ov[k] && k!="NM" && k!="MD" && k!="cs") bad[k]=1
            for(k in ov) if (!(k in iv)) bad[k]=1
            for(k in iv) if (!(k in ov)) bad[k]=1
          } END{ for(k in bad) printf "%s ", k }' "$T/${rt}.pair")"

  # NM/MD must agree with an independent implementation.
  samtools calmd -b "$out" "$REF" 2>/dev/null > "$T/${rt}.calmd.bam"
  tags_of() { samtools view "$1" | awk '{nm="";md=""; for(i=12;i<=NF;i++){if($i~/^NM:i:/)nm=$i; if($i~/^MD:Z:/)md=$i} print $1"_"$2"\t"nm"\t"md}' | sort; }
  paste <(tags_of "$out") <(tags_of "$T/${rt}.calmd.bam") > "$T/${rt}.calmd.cmp"
  chk "$rt NM agrees with samtools calmd" 0 \
      "$(awk -F'\t' '$1==$4 && $2!=$5 {n++} END{print n+0}' "$T/${rt}.calmd.cmp")"
  # Only compare MD where we emit one; MD is optional and absent from the HiFi input.
  chk "$rt MD agrees with samtools calmd" 0 \
      "$(awk -F'\t' '$1==$4 && $3!="" && $3!=$6 {n++} END{print n+0}' "$T/${rt}.calmd.cmp")"

  # Methylation calls must survive untouched.
  mm_of() { samtools view "$1" | awk '{m=""; for(i=12;i<=NF;i++) if($i~/^(MM|ML):/) m=m" "$i; print $1"_"$2"\t"m}' | sort; }
  chk "$rt MM/ML tags unchanged" same \
      "$([[ "$(mm_of "$BAM" | md5sum)" == "$(mm_of "$out" | md5sum)" ]] && echo same || echo differs)"
done

# ---------------------------------------------------------------------------
# Output completeness: reads on excluded contigs and unmapped reads.
#
# The bundled test files are chr11-only with no unmapped reads, so the pipeline
# visits every record and the gap is invisible. Synthesise a file that does have
# both, since without the passthrough pass the output BAM is a region-restricted
# subset rather than a drop-in replacement for the input.
# ---------------------------------------------------------------------------
HIFI=$ROOT/test_data/HG002_chr11_hifi_test.bam

mk_mixed() {  # $1=output bam; excluded contigs sort AFTER kept ones (normal GRCh38 layout)
  { samtools view --no-PG -H "$HIFI"
    samtools view "$HIFI"
    samtools view "$HIFI" | head -5 | awk 'BEGIN{OFS="\t"}{$1=$1"_M"; $3="chrM"; $4=100+NR*20; $5=60; print}'
    samtools view "$HIFI" | head -5 | awk 'BEGIN{OFS="\t"}{$1=$1"_D"; $3="chr1_KI270706v1_random"; $4=200+NR*20; $5=60; print}'
    samtools view "$HIFI" | head -3 | awk 'BEGIN{OFS="\t"}{$1=$1"_U"; $2=4; $3="*"; $4=0; $5=0; $6="*"; $7="*"; $8=0; $9=0; print}'
  } | samtools sort -o "$1" - 2>/dev/null
  samtools index "$1"
}

mk_reversed() {  # $1=output bam; chrM declared BEFORE chr11, so appending it cannot stay sorted
  { printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrM\tLN:16569\n@SQ\tSN:chr11\tLN:135086622\n'
    samtools view "$HIFI" | head -5 | awk 'BEGIN{OFS="\t"}{$1=$1"_M"; $3="chrM"; $4=100+NR*20; $5=60; print}'
    samtools view "$HIFI"
  } | samtools sort -o "$1" - 2>/dev/null
  samtools index "$1"
}

mk_mixed "$T/mixed.bam"
n_mixed=$(samtools view -c "$T/mixed.bam")

"$BIN" call --refine-aln "$REF" "$T/mixed.bam" -t 1 -b "$T/pass.bam" -o /dev/null >/dev/null 2>&1
chk "passthrough keeps every input record" "$n_mixed" "$(samtools view -c "$T/pass.bam")"
chk "passthrough output still coordinate-sorted" 0 "$(samtools view "$T/pass.bam" | awk '$3!="*"{if($3==p3 && $4<p4) n++; p3=$3;p4=$4} END{print n+0}')"
chk "passthrough output says SO:coordinate" "SO:coordinate" "$(samtools view --no-PG -H "$T/pass.bam" | grep -oE 'SO:[a-z]+' | head -1)"
rm -f "$T/pass.bam.bai"; samtools index "$T/pass.bam" 2>/dev/null
chk "passthrough output indexes" yes "$([[ -f $T/pass.bam.bai ]] && echo yes || echo no)"
chk "excluded-contig records are present" "10" \
    "$(( $(samtools view -c "$T/pass.bam" chrM) + $(samtools view -c "$T/pass.bam" chr1_KI270706v1_random) ))"
chk "unmapped records are present" 3 "$(samtools view -c -f 4 "$T/pass.bam")"
chk "passthrough records carry no HP/PS" 0 \
    "$(samtools view "$T/pass.bam" chrM chr1_KI270706v1_random | grep -cE 'HP:i:|PS:i:')"
chk "unmapped records carry no HP/PS" 0 "$(samtools view -f 4 "$T/pass.bam" | grep -cE 'HP:i:|PS:i:')"
chk "unmapped records come last" 0 "$(samtools view "$T/pass.bam" | awk '{if($3=="*") seen=1; else if(seen) n++} END{print n+0}')"
chk "qname+flag multiset preserved" same \
    "$([[ "$(samtools view "$T/mixed.bam" | awk '{print $1"_"$2}' | sort | md5sum)" == "$(samtools view "$T/pass.bam" | awk '{print $1"_"$2}' | sort | md5sum)" ]] && echo same || echo differs)"

# Opting out must drop exactly those records and stay sorted.
"$BIN" call --refine-aln --no-out-aln-passthrough "$REF" "$T/mixed.bam" -t 1 -b "$T/nopass.bam" -o /dev/null >/dev/null 2>&1
chk "--no-out-aln-passthrough omits them" "$((n_mixed - 13))" "$(samtools view -c "$T/nopass.bam")"

# A header that interleaves excluded and kept contigs cannot be both complete
# and sorted; keep the reads and say so.
mk_reversed "$T/rev.bam"
n_rev=$(samtools view -c "$T/rev.bam")
"$BIN" call --refine-aln "$REF" "$T/rev.bam" -t 1 -b "$T/rev.out.bam" -o /dev/null >/dev/null 2>&1
chk "interleaved header: records still all kept" "$n_rev" "$(samtools view -c "$T/rev.out.bam")"
chk "interleaved header: declares SO:unsorted" "SO:unsorted" "$(samtools view --no-PG -H "$T/rev.out.bam" | grep -oE 'SO:[a-z]+' | head -1)"

echo
if [[ $fail -eq 0 ]]; then echo "  BAM output: all $pass assertions passed"; exit 0
else echo "  BAM output: $fail of $((pass+fail)) assertions FAILED"; exit 1; fi
