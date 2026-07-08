#!/usr/bin/env bash
# Over-splitting test for a candidate population pair, via stepping-stone model
# evidence. Runs BA3 --ss on the K-population labeling and on the (K-1)-population
# labeling with popB merged into popA, and reports
#     log BF = logEv(K pops) - logEv(K-1 pops, popA+popB merged).
# BF < 0  => the split is not supported: popA, popB are over-split -> MERGE them.
# BF > 0  => popA, popB are genuinely distinct -> KEEP separate.
#
# The stepping-stone estimate has variance; the *sign* is stable but the magnitude
# is noisy for weak signals / small data. On real panels the BFs are large and
# unambiguous. If a result is close to 0, rerun with more iters and/or a second seed.
#
# Usage:
#   sim/ss_testmerge.sh <vcf> <meta> <popA> <popB> [iters] [seed]
#
# Env overrides:
#   BA3     path to the BA3 binary        (default: repo-root ./BA3)
#   GAMMA   "-G" fixes gamma=1 (default; the sensible choice here, since gamma is
#           unidentifiable for over-split pairs); set GAMMA="" to estimate it.
#   KEEP=1  keep the temp metadata/output files (printed path).
set -euo pipefail

if [ $# -lt 4 ]; then
	echo "usage: sim/ss_testmerge.sh <vcf> <meta> <popA> <popB> [iters=320000] [seed=7]" >&2
	exit 1
fi
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
vcf="$1"; meta="$2"; popA="$3"; popB="$4"
iters="${5:-320000}"; seed="${6:-7}"
BA3="${BA3:-$root/BA3}"
GAMMA="${GAMMA--G}"   # default: fix gamma=1

[ -x "$BA3" ] || { echo "error: BA3 not found/executable at $BA3 (build it, or set BA3=...)" >&2; exit 1; }
[ -f "$vcf" ]  || { echo "error: vcf not found: $vcf" >&2; exit 1; }
[ -f "$meta" ] || { echo "error: meta not found: $meta" >&2; exit 1; }
for p in "$popA" "$popB"; do
	awk -v p="$p" '$2==p{f=1} END{exit !f}' "$meta" || { echo "error: population '$p' not found in $meta" >&2; exit 1; }
done
K=$(awk 'NF>=2{print $2}' "$meta" | sort -u | wc -l | tr -d ' ')
if [ "$K" -lt 3 ]; then
	echo "error: merging leaves $((K-1)) population(s); BA3 needs >=2. The K=2 case" >&2
	echo "       (merge -> a single panmictic population) is not yet handled by this driver." >&2
	exit 1
fi

work="$(mktemp -d)"
cleanup() { [ "${KEEP:-0}" = 1 ] && echo "[keep] temp files in $work" || rm -rf "$work"; }
trap cleanup EXIT

# (K-1)-population metadata: popB -> popA
awk -v a="$popA" -v b="$popB" 'BEGIN{OFS="\t"} NF>=2{p=$2; if(p==b)p=a; print $1,p}' "$meta" > "$work/merged.txt"

ev() {  # <meta> -> stepping-stone log evidence
	"$BA3" -S $GAMMA -V "$vcf" -M "$1" -i "$iters" -b 0 -n 50 -s "$seed" -o "$work/out.txt" 2>/dev/null \
		| grep -oE 'log marginal likelihood = [-0-9.]+' | grep -oE '[-0-9.]+$'
}

echo "[ss] K=$K populations; testing merge of '$popB' into '$popA'"
echo "     (gamma: ${GAMMA:-estimated}, $iters iters/run, $seed seed, 16 rungs)"
eK=$(ev "$meta");            echo "  logEv($K pop)             = $eK"
eKm=$(ev "$work/merged.txt"); echo "  logEv($((K-1)) pop, merged)   = $eKm"
python3 - "$eK" "$eKm" "$popA" "$popB" <<'PY'
import sys
eK, eKm, a, b = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4]
bf = eK - eKm
verdict = (f"OVER-SPLIT -> merge {a} + {b}" if bf < 0
           else f"genuinely distinct -> keep {a} and {b} separate")
print(f"  log BF (split vs merged)  = {bf:+.1f}   ->  {verdict}")
PY
