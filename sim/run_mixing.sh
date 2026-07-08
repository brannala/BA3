#!/usr/bin/env bash
# Mixing benchmark: run BA3 in standard and --collapse modes on the same dataset
# with a dense (every-iteration) trace, then report integrated autocorrelation
# time (IAT) and effective sample size (ESS) of the migration-rate columns for
# each, via sim/ess.py. Quantifies the collapsed sampler's mixing advantage.
#
# Usage:
#   sim/run_mixing.sh <vcf> <meta> [iters] [burnin] [seed]
#
# Use a long chain: ESS/IAT estimates need many samples (n >> IAT) to be stable
# -- the defaults (100k iters, 20k burn-in) are a sensible floor; a few-thousand-
# sample run gives noisy, unreliable ESS and should not be compared.
#
# Env overrides:
#   BA3       path to the BA3 binary   (default: repo-root ./BA3)
#   PY        python with numpy        (default: autodetect, else build a venv)
#   KEEP=1    keep the large trace files (default: delete on exit)
set -euo pipefail

if [ $# -lt 2 ]; then
	echo "usage: sim/run_mixing.sh <vcf> <meta> [iters=100000] [burnin=20000] [seed=7]" >&2
	exit 1
fi

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
vcf="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
meta="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
iters="${3:-100000}"; burnin="${4:-20000}"; seed="${5:-7}"
BA3="${BA3:-$root/BA3}"

[ -x "$BA3" ] || { echo "error: BA3 binary not found/executable at $BA3 (build it, or set BA3=...)" >&2; exit 1; }

# Pick a python with numpy; otherwise build a cached venv next to this script.
if [ -z "${PY:-}" ]; then
	if python3 -c 'import numpy' 2>/dev/null; then
		PY=python3
	else
		# cache the venv outside the repo so it is not tracked or re-created each run
		venv="${XDG_CACHE_HOME:-$HOME/.cache}/ba3-mixvenv"
		[ -x "$venv/bin/python" ] || { echo "[setup] creating numpy venv at $venv" >&2; python3 -m venv "$venv" && "$venv/bin/pip" install -q numpy; }
		PY="$venv/bin/python"
	fi
fi

work="$(mktemp -d)"
cleanup() { [ "${KEEP:-0}" = 1 ] || rm -rf "$work"; }
trap cleanup EXIT

run() {  # run(mode_flag, out_trace_name)
	local flag="$1" name="$2"
	# BA3 writes BA3trace.txt to the CWD, so run inside the work dir.
	( cd "$work" && "$BA3" -V "$vcf" -M "$meta" $flag -t -i "$iters" -b "$burnin" -n 1 -s "$seed" -o "$work/${name}.out" >/dev/null 2>&1 )
	mv "$work/BA3trace.txt" "$work/${name}.trace"
}

echo "[run] standard  (i=$iters b=$burnin seed=$seed)  $(basename "$vcf")"
run ""   std
echo "[run] collapse  (--collapse)"
run "-c" col

echo
"$PY" "$here/ess.py" "standard" "$work/std.trace" "$burnin"
"$PY" "$here/ess.py" "collapse" "$work/col.trace" "$burnin"
[ "${KEEP:-0}" = 1 ] && echo "[keep] traces in $work"
echo "(ESS/iter is mixing per iteration; multiply by iters/sec for effective samples per second.)"
