#!/usr/bin/env python3
"""Convert BA3 format to VCF and metadata files."""

import sys
from collections import OrderedDict

def main():
    if len(sys.argv) != 4:
        print("Usage: ba3_to_vcf.py <input.ba3> <output.vcf> <output_meta.txt>")
        sys.exit(1)

    ba3_file = sys.argv[1]
    vcf_file = sys.argv[2]
    meta_file = sys.argv[3]

    # First pass: collect all individuals, populations, and loci
    indiv_to_pop = OrderedDict()
    loci = OrderedDict()

    print("First pass: collecting individuals, populations, and loci...")
    with open(ba3_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            indiv, pop, locus = parts[0], parts[1], parts[2]

            if indiv not in indiv_to_pop:
                indiv_to_pop[indiv] = pop

            if locus not in loci:
                loci[locus] = {}

    individuals = list(indiv_to_pop.keys())
    locus_names = list(loci.keys())

    print(f"Found {len(individuals)} individuals, {len(locus_names)} loci")

    # Initialize genotype storage: loci[locus][indiv] = (allele1, allele2)
    for locus in locus_names:
        loci[locus] = {ind: (None, None) for ind in individuals}

    # Second pass: read genotypes
    print("Second pass: reading genotypes...")
    with open(ba3_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            indiv, pop, locus, a1, a2 = parts[0], parts[1], parts[2], parts[3], parts[4]
            loci[locus][indiv] = (a1, a2)

    # Collect alleles per locus and determine REF/ALT
    print("Determining REF/ALT alleles per locus...")
    locus_alleles = {}
    for locus in locus_names:
        alleles = set()
        for indiv in individuals:
            a1, a2 = loci[locus][indiv]
            if a1 and a1 != '0':
                alleles.add(a1)
            if a2 and a2 != '0':
                alleles.add(a2)

        # Sort alleles alphabetically, first one is REF
        alleles = sorted(alleles)
        if len(alleles) == 0:
            alleles = ['N']  # All missing
        locus_alleles[locus] = alleles

    # Write metadata file
    print(f"Writing metadata to {meta_file}...")
    with open(meta_file, 'w') as f:
        for indiv, pop in indiv_to_pop.items():
            f.write(f"{indiv}\t{pop}\n")

    # Write VCF file
    print(f"Writing VCF to {vcf_file}...")
    with open(vcf_file, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=1000000000>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        # Column header
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        header.extend(individuals)
        f.write("\t".join(header) + "\n")

        # Data lines
        pos = 0
        for locus in locus_names:
            pos += 1
            alleles = locus_alleles[locus]
            ref = alleles[0]
            alt = ",".join(alleles[1:]) if len(alleles) > 1 else "."

            # Create allele-to-index mapping
            allele_idx = {a: str(i) for i, a in enumerate(alleles)}

            # Build genotype strings
            gts = []
            for indiv in individuals:
                a1, a2 = loci[locus][indiv]

                # Handle missing data
                if a1 == '0' or a1 is None:
                    gt1 = "."
                else:
                    gt1 = allele_idx.get(a1, ".")

                if a2 == '0' or a2 is None:
                    gt2 = "."
                else:
                    gt2 = allele_idx.get(a2, ".")

                gts.append(f"{gt1}/{gt2}")

            # Parse locus name for chromosome info (e.g., RAD_0-8 -> chr0, pos 8)
            chrom = "chr1"

            row = [chrom, str(pos), locus, ref, alt, ".", "PASS", ".", "GT"]
            row.extend(gts)
            f.write("\t".join(row) + "\n")

    print("Done!")

if __name__ == "__main__":
    main()
