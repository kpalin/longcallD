#!/usr/bin/env bash
# Parameter sweep for paired calling: run the real caller once per config and
# score every result with som.py, so the numbers are comparable to the headline
# benchmark rather than to a reimplementation of the classifier.
#
# That distinction is not pedantic. An earlier attempt replayed the four-slot
# decision offline from the VCF's own HDP/HAD evidence, which is far cheaper,
# and it was useless: it reported 1863 TP / 222 FP where som.py measured
# 2108 / 76, because it omitted the germline and LOH branches and the pooled
# fallback. Only the real binary knows what the real binary does.
#
# Scoring is restricted to the sweep region with -l. Leaving it at chr1 makes
# the whole chromosome's truth set the denominator while the query covers only
# the window, which silently pins FN near 3400 and makes recall and F1
# meaningless while leaving TP, FP and precision looking perfectly reasonable.
#
# Usage: eval/sweep_params.sh -c CONFIGS [-r REGION] [-t THREADS] [-o OUTDIR]
#
#   -c  config file: one "name|args" per line, '#' comments and blanks ignored.
#       An empty args field is the defaults baseline.
#
# Each config's VCF is kept, so a re-run re-scores without re-calling: handy
# because scoring is seconds and calling is minutes.
set -uo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DATA=${DATA:-/home/ubuntu/castle_tn}
REF=${REF:-$DATA/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna}
TUMOR=${TUMOR:-$DATA/HCC1395_pacbio.tumor.chr1.bam}
NORMAL=${NORMAL:-$DATA/HCC1395_pacbio.normal.chr1.bam}
BIN=${BIN:-$ROOT/bin/longcallD}
MEMCAP_KB=${MEMCAP_KB:-5000000}
REGION=chr1:1-25000000; THREADS=2; OUT=${OUT:-/tmp/claude-1000/sweep}; CONFIGS=""
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
while getopts "c:r:t:o:h" o; do case $o in
  c) CONFIGS=$OPTARG ;; r) REGION=$OPTARG ;; t) THREADS=$OPTARG ;; o) OUT=$OPTARG ;;
  *) sed -n '2,20p' "$0"; exit 0 ;;
esac; done
[[ -n $CONFIGS && -s $CONFIGS ]] || die "need -c CONFIGS (see -h)"
mkdir -p "$OUT"
printf "%-24s %5s %5s %5s %6s %7s %7s %7s\n" CONFIG TP FP FN NCALL RECALL PREC F1
while IFS='|' read -r name args; do
  [[ -z ${name// /} || ${name:0:1} == '#' ]] && continue
  name=${name// /}; vcf="$OUT/$name.vcf"
  if [[ ! -s $vcf ]]; then
    # shellcheck disable=SC2086
    ( ulimit -v "$MEMCAP_KB"
      "$BIN" paired -t "$THREADS" $args "$REF" "$TUMOR" "$NORMAL" "$REGION" -o "$vcf.tmp" ) \
      >/dev/null 2>"$OUT/$name.log" \
      && mv "$vcf.tmp" "$vcf" || { printf "%-24s CALL FAILED (see %s)\n" "$name" "$OUT/$name.log"; rm -f "$vcf.tmp"; continue; }
  fi
  # The scored query is the PASS somatic set, exactly as run_benchmark.sh builds
  # it -- germline calls also PASS and must come out.
  bcftools view -f PASS -i 'INFO/SOMATIC=1' -Oz -o "$OUT/$name.q.vcf.gz" "$vcf" 2>/dev/null
  tabix -f -p vcf "$OUT/$name.q.vcf.gz" 2>/dev/null
  ncall=$(bcftools view -H "$OUT/$name.q.vcf.gz" | wc -l)
  "$ROOT/eval/sompy.sh" -q "$OUT/$name.q.vcf.gz" -o "$OUT/$name.sompy" -l "$REGION" >"$OUT/$name.score" 2>&1 \
    || { printf "%-24s SCORE FAILED (see %s)\n" "$name" "$OUT/$name.score"; continue; }
  # The env prefix must carry S as well: a VAR=value written after the command
  # is an argument to python, not an assignment, so S was never in os.environ.
  NAME=$name NCALL=$ncall S="$OUT/$name.sompy.stats.csv" python3 -I -c "
import csv, os
rows=[r for r in csv.DictReader(open(os.environ['S'])) if r['type']=='records']
if not rows: print('%-24s no records row'%os.environ['NAME']); raise SystemExit
r=rows[0]; tp,fp,fn=int(r['tp']),int(r['fp']),int(r['fn'])
rec=tp/(tp+fn) if tp+fn else 0; pre=tp/(tp+fp) if tp+fp else 0
f1=2*rec*pre/(rec+pre) if rec+pre else 0
print('%-24s %5d %5d %5d %6s %7.4f %7.4f %7.4f'%(os.environ['NAME'],tp,fp,fn,os.environ['NCALL'],rec,pre,f1))
"
done < "$CONFIGS"
