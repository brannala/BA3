#!/usr/bin/env bash
# Regression test: the standard (non --collapse) sampler must be byte-identical
# to a reference BA3 build for the same seed. Every output file (results, trace,
# per-individual, allele-frequency) and stdout is compared; only the wall-clock
# "Elapsed time" line and the progress bar's ETA are excluded.
#
# usage: tests/regress_default.sh <reference BA3> <new BA3> [iters=20000]
#   e.g. build master into /tmp/ref and run
#        tests/regress_default.sh /tmp/ref/BA3 ./BA3
set -uo pipefail

REF="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
NEW="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
ITERS="${3:-20000}"
BURN=$((ITERS / 10))
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT

# name | input args (native file, or -V/-M for VCF)
CASES=(
	"2pop|$root/examples/2pop.txt"
	"3pop|$root/examples/3pop.txt"
	"100loci_mig|$root/examples/100loci_with_migration.txt"
	"high_info|$root/examples/high_info_no_migration.txt"
	"2pop_missing|$root/tests/data/2pop_missing.txt"
	"100loci_missing|$root/tests/data/100loci_missing.txt"
	"vcf_litt|-V $root/examples/allpopslitt.vcf -M $root/examples/allpopslitt_meta.txt"
)
# extra option sets exercised on every case
OPTS=(
	""
	"-t -g -F freq.txt"
	"-N -v"
)

run() {  # run <bin> <dir> <opts> <input args>
	mkdir -p "$2"
	( cd "$2" && "$1" $3 -s 7 -i "$ITERS" -b "$BURN" -n 50 -o out.txt $4 \
		2>&1 | tr '\r' '\n' | sed -E 's/ETA: .*//' | grep -v 'Elapsed time' > stdout.txt )
}

fail=0; n=0
for c in "${CASES[@]}"; do
	name="${c%%|*}"; input="${c#*|}"
	for k in "${!OPTS[@]}"; do
		o="${OPTS[$k]}"; d="$work/$name.$k"
		run "$REF" "$d/ref" "$o" "$input"
		run "$NEW" "$d/new" "$o" "$input"
		n=$((n + 1))
		if diff -r -q "$d/ref" "$d/new" > "$d/diff.txt"; then
			printf "  PASS  %-16s opts=[%s]\n" "$name" "$o"
		else
			printf "  FAIL  %-16s opts=[%s]\n" "$name" "$o"; sed 's/^/        /' "$d/diff.txt"
			fail=$((fail + 1))
		fi
	done
done
echo "regress_default: $((n - fail))/$n identical"
[ "$fail" -eq 0 ]
