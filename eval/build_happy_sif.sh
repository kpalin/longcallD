#!/usr/bin/env bash
#
# Build an apptainer SIF of pkrusche/hap.py -- the exact container the
# DeepSomatic paper used for som.py (hap.py v0.3.9, pushed 2017-06-03).
#
# WHY THIS IS NOT JUST `apptainer pull`:
#   $ apptainer pull happy.sif docker://pkrusche/hap.py
#   FATAL: unsupported MediaType:
#          "application/vnd.docker.distribution.manifest.v1+prettyjws"
#
# Every tag of that image (latest, v0.3.6 .. v0.3.9) is stored as a Docker
# manifest *schema v1*, and Docker Hub serves schema v1 no matter what the
# client puts in its Accept header -- verified against registry-1.docker.io.
# apptainer's go-containerregistry backend dropped schema v1 support
# (google/go-containerregistry#377), so no modern OCI client can pull it.
#
# skopeo still reads schema v1 and rewrites it as OCI on the way out, so the
# conversion is a two-step: skopeo copy -> apptainer build. The image content
# is unchanged; only the manifest format is.
#
set -euo pipefail

DATA=${DATA:-/home/ubuntu/castle_tn}
IMAGE=${IMAGE:-docker://pkrusche/hap.py:latest}
SIF=${SIF:-$DATA/happy.sif}
OCI=${OCI:-$DATA/happy.oci}

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

command -v skopeo    >/dev/null || die "skopeo not found. sudo apt-get install -y skopeo"
command -v apptainer >/dev/null || die "apptainer not found"

# Rootless apptainer needs unprivileged user namespaces.
if [[ -r /proc/sys/user/max_user_namespaces ]]; then
  [[ $(< /proc/sys/user/max_user_namespaces) -gt 0 ]] \
    || die "unprivileged user namespaces are disabled; rootless apptainer will not work"
fi

mkdir -p "$DATA"

if [[ -s $SIF ]]; then
  log "SIF already present: $SIF"
else
  if [[ ! -s $OCI ]]; then
    log "skopeo copy $IMAGE -> $OCI (schema v1 -> OCI, ~970 MB)"
    skopeo copy "$IMAGE" "oci-archive:$OCI:happy"
  else
    log "OCI archive already present: $OCI"
  fi
  log "apptainer build $SIF (~525 MB)"
  APPTAINER_CACHEDIR="$DATA/.apptainer_cache" apptainer build --force "$SIF" "oci-archive:$OCI"
fi

log "verifying som.py runs"
apptainer exec "$SIF" /opt/hap.py/bin/som.py --help >/dev/null 2>&1 \
  || die "som.py did not run inside $SIF"
log "ok: $SIF"
