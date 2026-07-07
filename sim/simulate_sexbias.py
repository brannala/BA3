#!/usr/bin/env python3
"""
simulate_sexbias.py

Simple forward simulator producing ground-truth data for testing BA3's
sex-biased dispersal model (see doc/sex_biased_dispersal_model.md).

It generates diploid autosomal and X-linked biallelic SNP genotypes for a set of
populations connected by known, sex-specific migration, then writes:

  <out>.vcf        VCF (chr1 = autosomal loci, chrX = X-linked loci; males are
                   haploid/hemizygous on chrX)
  <out>_meta.txt   metadata, tab-separated: sampleID <tab> pop <tab> sex(F/M)
  <out>_truth.json ground-truth parameters and per-individual true ancestry

Model (matches the design doc, migrantAge in {0,1,2}):
  * Each individual is a non-migrant (age 0), a first-generation migrant from a
    source population (age 1), or a second-generation migrant / F1 hybrid whose
    single migrant parent came from a source (age 2).
  * The first-generation migrant is female with probability phi (the global
    dispersal female fraction). A first-gen migrant (age 1) *is* the sampled
    individual, so its own sex equals the migrant sex. A second-gen migrant
    (age 2) has its own sex drawn 1/2, and its migrant parent's sex ~ Bernoulli(phi).
  * Genotypes are drawn under HWE (F = 0) from per-population allele frequencies.
    X inheritance follows the copy-type rules of the design doc:
       - age 1 : both / single X copy from the source population
       - age 2 male   : X all-source if migrant parent is the mother, else
                        all-native (a male's single X comes from his mother)
       - age 2 female : one X from source, one from native (regardless of which
                        parent migrated)

Population differentiation is generated with a Balding-Nichols draw controlled
by Fst. Everything is pure Python stdlib (no numpy required).

Example:
    python3 sim/simulate_sexbias.py --npop 2 --nind 60 --phi 0.8 \
        --m1 0.10 --m2 0.05 --nauto 300 --nx 150 --fst 0.1 --seed 42 --out sim_run
"""

import argparse
import json
import random
import sys


def draw_pop_freqs(ancestral, fst, npop, rng):
    """Balding-Nichols per-population ALT allele frequency for one locus.

    Returns a list of length npop. With fst == 0 there is no differentiation
    (every population equals the ancestral frequency)."""
    if fst <= 0.0:
        return [ancestral] * npop
    a = ancestral * (1.0 - fst) / fst
    b = (1.0 - ancestral) * (1.0 - fst) / fst
    return [rng.betavariate(a, b) for _ in range(npop)]


def build_freqs(nloci, fst, npop, rng, fmin):
    """For each locus draw an ancestral ALT frequency then per-population freqs.
    Returns freqs[locus][pop] = ALT allele frequency."""
    freqs = []
    for _ in range(nloci):
        ancestral = rng.uniform(fmin, 1.0 - fmin)
        freqs.append(draw_pop_freqs(ancestral, fst, npop, rng))
    return freqs


def copy(freq, rng):
    """Draw one gene copy: 1 (ALT) with prob freq, else 0 (REF)."""
    return 1 if rng.random() < freq else 0


def assign_ancestry(s, npop, m1, m2, rng):
    """Draw (age, source) for one individual in population s.

    age 0 -> (0, None); age 1/2 -> (age, source_pop)."""
    others = [b for b in range(npop) if b != s]
    # Categorical over: native, (age1, b) for each b, (age2, b) for each b.
    labels = [(0, None)]
    weights = [1.0 - (m1 + m2) * len(others)]
    for b in others:
        labels.append((1, b))
        weights.append(m1)
    for b in others:
        labels.append((2, b))
        weights.append(m2)
    r = rng.random()
    cumulative = 0.0
    for (label, w) in zip(labels, weights):
        cumulative += w
        if r < cumulative:
            return label
    return labels[-1]  # numerical safety


def simulate(args):
    rng = random.Random(args.seed)

    npop = args.npop
    total_mig_rate = (args.m1 + args.m2) * (npop - 1)
    if total_mig_rate >= 1.0:
        sys.exit("Error: (m1 + m2) * (npop - 1) = %.3f must be < 1." % total_mig_rate)

    # Per-population allele frequencies (ALT), separately for autosomal and X.
    auto_freq = build_freqs(args.nauto, args.fst, npop, rng, args.fmin)
    x_freq = build_freqs(args.nx, args.fst, npop, rng, args.fmin)

    individuals = []      # metadata / truth per individual
    auto_geno = []        # auto_geno[i] = list of (a0, a1) over autosomal loci
    x_geno = []           # x_geno[i]    = list of (a0, a1) (female) or (a0,) (male)

    for s in range(npop):
        for k in range(args.nind):
            age, source = assign_ancestry(s, npop, args.m1, args.m2, rng)

            # Sex assignment and the migrant's sex.
            if age == 0:
                sex = "F" if rng.random() < 0.5 else "M"
                migrant_sex = None
            elif age == 1:
                # The sampled individual IS the migrant.
                migrant_sex = "F" if rng.random() < args.phi else "M"
                sex = migrant_sex
            else:  # age == 2
                sex = "F" if rng.random() < 0.5 else "M"
                migrant_sex = "F" if rng.random() < args.phi else "M"

            sid = "pop%d_%d" % (s, k)
            individuals.append({
                "id": sid, "pop": "pop%d" % s, "pop_index": s, "sex": sex,
                "age": age, "source": None if source is None else "pop%d" % source,
                "source_index": source, "migrant_sex": migrant_sex,
            })

            # ---- autosomal genotypes ----
            g_auto = []
            for l in range(args.nauto):
                fs = auto_freq[l][s]
                if age == 0:
                    a0, a1 = copy(fs, rng), copy(fs, rng)
                elif age == 1:
                    fb = auto_freq[l][source]
                    a0, a1 = copy(fb, rng), copy(fb, rng)
                else:  # age 2 hybrid: one source copy, one native copy
                    fb = auto_freq[l][source]
                    a0, a1 = copy(fb, rng), copy(fs, rng)
                g_auto.append((a0, a1))
            auto_geno.append(g_auto)

            # ---- X-linked genotypes ----
            g_x = []
            for l in range(args.nx):
                fs = x_freq[l][s]
                fb = x_freq[l][source] if source is not None else None
                if sex == "F":
                    if age == 0:
                        a0, a1 = copy(fs, rng), copy(fs, rng)
                    elif age == 1:
                        a0, a1 = copy(fb, rng), copy(fb, rng)
                    else:  # age 2 female: one source X, one native X
                        a0, a1 = copy(fb, rng), copy(fs, rng)
                    g_x.append((a0, a1))
                else:  # male: single (hemizygous) X copy, from his mother
                    if age == 0:
                        a0 = copy(fs, rng)
                    elif age == 1:
                        a0 = copy(fb, rng)
                    else:  # age 2 male: source X iff migrant parent is the mother
                        a0 = copy(fb, rng) if migrant_sex == "F" else copy(fs, rng)
                    g_x.append((a0,))
            x_geno.append(g_x)

    return individuals, auto_geno, x_geno, auto_freq, x_freq


def write_vcf(path, individuals, auto_geno, x_geno, nauto, nx):
    sample_ids = [ind["id"] for ind in individuals]
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=simulate_sexbias.py\n")
        fh.write("##contig=<ID=chr1,length=%d>\n" % max(nauto, 1))
        fh.write("##contig=<ID=chrX,length=%d>\n" % max(nx, 1))
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        fh.write("\t".join(header + sample_ids) + "\n")

        # Autosomal rows (diploid for everyone).
        for l in range(nauto):
            row = ["chr1", str(l + 1), "auto_%d" % (l + 1), "A", "C", ".", "PASS", ".", "GT"]
            for i in range(len(individuals)):
                a0, a1 = auto_geno[i][l]
                row.append("%d/%d" % (a0, a1))
            fh.write("\t".join(row) + "\n")

        # X rows (females diploid, males haploid/hemizygous).
        for l in range(nx):
            row = ["chrX", str(l + 1), "x_%d" % (l + 1), "A", "C", ".", "PASS", ".", "GT"]
            for i in range(len(individuals)):
                gt = x_geno[i][l]
                if len(gt) == 2:
                    row.append("%d/%d" % (gt[0], gt[1]))
                else:
                    row.append("%d" % gt[0])
            fh.write("\t".join(row) + "\n")


def write_metadata(path, individuals):
    with open(path, "w") as fh:
        for ind in individuals:
            fh.write("%s\t%s\t%s\n" % (ind["id"], ind["pop"], ind["sex"]))


def summarize(individuals):
    """Empirical checks: realized female fraction among first-gen migrants and
    per-category counts."""
    n_fg_female = n_fg_male = 0
    counts = {"age0": 0, "age1": 0, "age2": 0}
    for ind in individuals:
        counts["age%d" % ind["age"]] += 1
        if ind["migrant_sex"] == "F":
            n_fg_female += 1
        elif ind["migrant_sex"] == "M":
            n_fg_male += 1
    total_mig = n_fg_female + n_fg_male
    emp_phi = (n_fg_female / total_mig) if total_mig else float("nan")
    return {
        "counts_by_age": counts,
        "n_migrant_lineages": total_mig,
        "n_first_gen_migrant_female": n_fg_female,
        "n_first_gen_migrant_male": n_fg_male,
        "empirical_phi": emp_phi,
    }


def main():
    p = argparse.ArgumentParser(description="Simulate sex-biased dispersal data for BA3 testing.")
    p.add_argument("--npop", type=int, default=2, help="number of populations")
    p.add_argument("--nind", type=int, default=60, help="individuals per population")
    p.add_argument("--phi", type=float, default=0.8,
                   help="female fraction of migrants (0..1); 0.5 = no sex bias")
    p.add_argument("--m1", type=float, default=0.10,
                   help="per-source first-generation migrant rate (age 1)")
    p.add_argument("--m2", type=float, default=0.05,
                   help="per-source second-generation migrant rate (age 2)")
    p.add_argument("--nauto", type=int, default=300, help="number of autosomal SNP loci")
    p.add_argument("--nx", type=int, default=150, help="number of X-linked SNP loci")
    p.add_argument("--fst", type=float, default=0.1,
                   help="Balding-Nichols Fst controlling population differentiation")
    p.add_argument("--fmin", type=float, default=0.05,
                   help="minimum ancestral ALT frequency (avoids monomorphic loci)")
    p.add_argument("--seed", type=int, default=42, help="RNG seed")
    p.add_argument("--out", type=str, default="sim_run", help="output file prefix")
    args = p.parse_args()

    if not (0.0 <= args.phi <= 1.0):
        sys.exit("Error: --phi must be in [0, 1].")

    individuals, auto_geno, x_geno, auto_freq, x_freq = simulate(args)

    vcf_path = args.out + ".vcf"
    meta_path = args.out + "_meta.txt"
    truth_path = args.out + "_truth.json"

    write_vcf(vcf_path, individuals, auto_geno, x_geno, args.nauto, args.nx)
    write_metadata(meta_path, individuals)

    summary = summarize(individuals)
    truth = {
        "params": {
            "npop": args.npop, "nind_per_pop": args.nind, "phi": args.phi,
            "m1": args.m1, "m2": args.m2, "nauto": args.nauto, "nx": args.nx,
            "fst": args.fst, "fmin": args.fmin, "seed": args.seed,
        },
        "summary": summary,
        "individuals": individuals,
        "auto_freq": auto_freq,
        "x_freq": x_freq,
    }
    with open(truth_path, "w") as fh:
        json.dump(truth, fh, indent=1)

    print("Wrote:")
    print("  %s   (%d autosomal + %d X loci, %d individuals)"
          % (vcf_path, args.nauto, args.nx, len(individuals)))
    print("  %s" % meta_path)
    print("  %s" % truth_path)
    print("Ground truth: phi = %.3f  (empirical among %d migrant lineages: %.3f)"
          % (args.phi, summary["n_migrant_lineages"], summary["empirical_phi"]))
    print("Ancestry counts: %s" % summary["counts_by_age"])


if __name__ == "__main__":
    main()
