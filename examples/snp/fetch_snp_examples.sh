#!/usr/bin/env bash
# Download the large SNP example datasets from Dryad (CC0) and prepare BA3 inputs
# (VCF + INDIV/POPLN metadata). Idempotent: skips downloads already present.
# Run from examples/snp/.  See README.md for provenance and citations.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

# Dryad public file-stream downloads (numeric file ids from the datasets' API).
DACE_URL="https://datadryad.org/stash/downloads/file_stream/394661"    # doi:10.5061/dryad.51c59zw62
PALM_URL="https://datadryad.org/stash/downloads/file_stream/3107434"   # doi:10.5061/dryad.pc866t1z6

fetch() {  # url outfile
	if [ -s "$2" ]; then echo "have $(basename "$2")"; return; fi
	echo "downloading $(basename "$2") ..."
	curl -fSL "$1" -o "$2" || { echo "  download failed; grab it in a browser from the DOI in README.md and re-run" >&2; exit 1; }
}

# ---- dace: 145 indiv, 13 pops, 14,355 SNPs; pop = letters in the sample name ----
mkdir -p dace
fetch "$DACE_URL" "$tmp/dace.tar.gz"
tar -xzf "$tmp/dace.tar.gz" -C "$tmp"
cp "$(find "$tmp" -name spd_dv.filt2.vcf | head -1)" dace/spd_dv.filt2.vcf
grep -m1 '^#CHROM' dace/spd_dv.filt2.vcf | cut -f10- | tr '\t' '\n' \
	| awk '{n=$0; p=$0; sub(/^[0-9]+/,"",p); sub(/[0-9]+$/,"",p); print n"\t"toupper(p)}' \
	> dace/spd_dv_meta.txt
echo "  dace: $(grep -vc '^#' dace/spd_dv.filt2.vcf) loci, $(wc -l < dace/spd_dv_meta.txt) indiv"

# ---- palm: 4 species, *_miss50.vcf; pop = sample-name prefix before first '_' ----
mkdir -p palm
fetch "$PALM_URL" "$tmp/palm_outer.zip"
unzip -oq "$tmp/palm_outer.zip" -d "$tmp/palm"
inner="$(find "$tmp/palm" -name '*.zip' | head -1)"
[ -n "$inner" ] && unzip -oq "$inner" -d "$tmp/palm/x" || cp -r "$tmp/palm" "$tmp/palm/x"
for sp in Bm BN Cm HC; do
	v="$(find "$tmp/palm/x" -name "${sp}_miss50.vcf" | head -1)"
	[ -n "$v" ] || { echo "  palm $sp: VCF not found" >&2; continue; }
	cp "$v" "palm/${sp}_miss50.vcf"
	grep -m1 '^#CHROM' "$v" | cut -f10- | tr '\t' '\n' | awk -F_ '{print $0"\t"$1}' > "palm/${sp}_meta.txt"
	echo "  palm $sp: $(grep -vc '^#' "palm/${sp}_miss50.vcf") loci, $(wc -l < "palm/${sp}_meta.txt") indiv"
done

echo "done. Example: ../../BA3 -c -V dace/spd_dv.filt2.vcf -M dace/spd_dv_meta.txt -i 5000000 -b 1000000 -n 2000 -o dace_col.txt"
