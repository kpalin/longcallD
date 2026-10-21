#!/usr/bin/env bash
#
# The per-PR C testing gate: build longcallD under each sanitizer and run the
# bundled test data through it, failing on any finding.
#
# Deliberately uses the bundled 2 Mb test_data/ and never chr1 -- ASan costs
# roughly 2x time and 3x RSS, so a sanitized chr1 run at 40x+55x would OOM on
# a small box long before it found anything.
#
# Two environment quirks this encodes, both found the hard way:
#
#   * TSan aborts with "unexpected memory mapping" under this kernel's ASLR
#     entropy, so TSan runs go through `setarch -R` to disable randomisation.
#   * klib's work stealing has a deliberate benign race; see eval/tsan.supp.
#
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
OUT=${OUT:-$(mktemp -d)}
mkdir -p "$OUT"
JOBS=${JOBS:-4}
REF=test_data/chr11_2M.fa
fails=0

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
ok()  { printf '  PASS  %s\n' "$1"; }
bad() { printf '  FAIL  %s\n' "$1"; fails=$((fails+1)); }

# A Makefile that does not understand san= builds a perfectly ordinary binary,
# and then every check below passes without testing anything. Refuse to report
# a pass we have not earned.
assert_instrumented() {   # $1=label $2=symbol-regex
  local n
  n=$(nm bin/longcallD 2>/dev/null | grep -cE "$2")
  if [[ ${n:-0} -lt 1 ]]; then
    bad "$1 -- bin/longcallD has no $2 symbols; does this Makefile support san= ?"
    return 1
  fi
  ok "$1 instrumentation present ($n symbols)"
  return 0
}

run_matrix() {   # $1=label $2=findings-pattern $3.. = launcher prefix
  local label=$1 pattern=$2; shift 2
  local rt t refine flag ec logf name
  for rt in hifi ont; do
    for t in 1 4; do
      for refine in "" "--refine-aln"; do
        flag=""; [[ $rt == ont ]] && flag="--ont"
        name="$label ${rt} -t${t}${refine:+ ${refine}}"
        logf="$OUT/${label}_${rt}_t${t}${refine:+_refine}.log"
        if [[ -n $refine ]]; then
          "$@" ./bin/longcallD call $flag --refine-aln "$REF" \
              "test_data/HG002_chr11_${rt}_test.bam" -t "$t" \
              -b "$OUT/${label}_${rt}_t${t}.bam" \
              -o "$OUT/${label}_${rt}_t${t}.vcf" >"$logf" 2>&1
        else
          "$@" ./bin/longcallD call $flag "$REF" \
              "test_data/HG002_chr11_${rt}_test.bam" -t "$t" \
              -o "$OUT/${label}_${rt}_t${t}.vcf" >"$logf" 2>&1
        fi
        ec=$?
        if grep -qE "$pattern" "$logf"; then
          bad "$name -- $(grep -oE "$pattern" "$logf" | head -1)"
        elif [[ $ec -ne 0 ]]; then
          bad "$name -- exit $ec"
        else
          ok "$name"
        fi
      done
    done
  done
}

log "=== ASan + UBSan (make san=1) ==="
make clean >/dev/null 2>&1
if ! make san=1 -j"$JOBS" >"$OUT/build_asan.log" 2>&1; then
  bad "san=1 build"; tail -3 "$OUT/build_asan.log"
else
  assert_instrumented asan '__asan|__ubsan'
  export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
  export ASAN_OPTIONS=detect_leaks=1
  run_matrix asan 'ERROR: (Address|Leak|Undefined)Sanitizer|runtime error|FATAL'
fi

log "=== TSan (make san=thread) ==="
make clean >/dev/null 2>&1
if ! make san=thread -j"$JOBS" >"$OUT/build_tsan.log" 2>&1; then
  bad "san=thread build"; tail -3 "$OUT/build_tsan.log"
else
  assert_instrumented tsan '__tsan'
  unset UBSAN_OPTIONS ASAN_OPTIONS
  export TSAN_OPTIONS="halt_on_error=0:suppressions=$ROOT/eval/tsan.supp"
  run_matrix tsan 'WARNING: ThreadSanitizer|FATAL' setarch "$(uname -m)" -R
fi

# Leave a normal build behind. Otherwise bin/longcallD is still the TSan binary,
# which exits 66 on any reported race, and the next unrelated test run fails in
# a way that looks like a code regression.
log "=== restoring a normal build ==="
make clean >/dev/null 2>&1
if make -j"$JOBS" >"$OUT/build_normal.log" 2>&1; then
  log "normal build restored"
else
  bad "could not restore the normal build"; tail -3 "$OUT/build_normal.log"
fi

echo
if [[ $fails -eq 0 ]]; then
  log "ALL SANITIZER CHECKS PASSED (logs in $OUT)"; exit 0
else
  log "$fails SANITIZER CHECK(S) FAILED (logs in $OUT)"; exit 1
fi
