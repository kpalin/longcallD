#!/usr/bin/env bash
#
# E0 -- fetch the HCC1395/HCC1395BL PacBio HiFi chr1 fixtures and the SEQC2
# somatic truth set used by the DeepSomatic paper (Park et al., Nat Biotechnol
# 2025, doi 10.1038/s41587-025-02839-x).
#
# These are the exact chr1 subsets the DeepSomatic PacBio case study evaluates,
# already aligned with pbmm2 to GRCh38 (chr-prefixed) and carrying no PacBio
# kinetics tags (only 5mC MM/ML), so nothing needs realigning or stripping.
#
# Resumable: curl -C - continues partial transfers, and any file already at its
# expected byte size is skipped. Safe to re-run after an interruption.
#
set -euo pipefail

DATA=${DATA:-/home/ubuntu/castle_tn}

GCS=https://storage.googleapis.com/deepvariant/deepsomatic-case-studies
CS=$GCS/deepsomatic-chr1-case-studies
TR=$GCS/SEQC2-S1395-truth

# url <TAB> expected_size   (sizes verified 2026-10-21; a mismatch is fatal,
# because a silently truncated BAM looks like a low-coverage sample)
read -r -d '' SMALL <<EOF || true
$TR/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz	5108710
$TR/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz.tbi	367972
$TR/High-Confidence_Regions_v1.2.bed	16879606
$CS/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna	253105702
$CS/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna.fai	23
EOF

read -r -d '' BIG <<EOF || true
$CS/HCC1395_pacbio.normal.chr1.bam	5136335408
$CS/HCC1395_pacbio.normal.chr1.bam.bai	1775368
$CS/HCC1395_pacbio.tumor.chr1.bam	6711493364
$CS/HCC1395_pacbio.tumor.chr1.bam.bai	3385152
EOF

usage() {
  cat <<EOF
Usage: $0 [-s] [-b] [-d DATA_DIR]

  -s   fetch only the small files (truth VCF, HC BED, chr1 reference) ~275 MB
  -b   fetch only the big files (tumor + normal chr1 BAMs)            ~11.9 GB
  -d   data directory (default: $DATA)

With neither -s nor -b, fetches everything (~12.1 GB).
EOF
  exit 1
}

do_small=0 do_big=0
while getopts "sbd:h" o; do
  case $o in
    s) do_small=1 ;;
    b) do_big=1 ;;
    d) DATA=$OPTARG ;;
    *) usage ;;
  esac
done
if [[ $do_small -eq 0 && $do_big -eq 0 ]]; then do_small=1; do_big=1; fi

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

command -v curl >/dev/null || die "curl not found"
mkdir -p "$DATA"

fetch() {
  local url=$1 want=$2 out="$DATA/${1##*/}"
  if [[ -f $out ]]; then
    local have
    have=$(stat -c %s "$out")
    if [[ $have -eq $want ]]; then
      log "ok (cached): ${url##*/}"
      return 0
    fi
    log "resuming ${url##*/} ($have / $want bytes)"
  else
    log "fetching ${url##*/} ($want bytes)"
  fi
  curl -fL --retry 5 --retry-delay 10 -C - -o "$out" "$url"
  local got
  got=$(stat -c %s "$out")
  [[ $got -eq $want ]] || die "${url##*/}: got $got bytes, expected $want"
  log "ok: ${url##*/}"
}

[[ $do_small -eq 1 ]] && while IFS=$'\t' read -r url want; do
  [[ -n ${url:-} ]] && fetch "$url" "$want"
done <<< "$SMALL"

[[ $do_big -eq 1 ]] && while IFS=$'\t' read -r url want; do
  [[ -n ${url:-} ]] && fetch "$url" "$want"
done <<< "$BIG"

log "done. data in $DATA"
