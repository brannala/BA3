#!/usr/bin/env python3
"""Validation gate for the collapsed (integrated-allele-frequency) sampler.

Phase 1: runs `BA3 --collapse` on a VCF+metadata dataset, reads the reported
"collapsed init logL" (the Dirichlet-multinomial marginal of the initial
all-native count tables), and checks it against an independent computation here.

Usage:
  python3 sim/check_collapsed.py ./BA3 data.vcf data_meta.txt
"""
import sys, math, subprocess, collections, re

ALPHA = 1.0


def independent_dirmult(vcf, meta):
    pop = {}
    for l in open(meta):
        p = l.split()
        if len(p) >= 2:
            pop[p[0]] = p[1]
    counts = collections.defaultdict(lambda: [0, 0])   # (pop, locus) -> [nref, nalt]; biallelic
    samples = None
    li = 0
    for line in open(vcf):
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
            continue
        if line.startswith('#'):
            continue
        gts = line.rstrip('\n').split('\t')[9:]
        for s, g in zip(samples, gts):
            for a in g.split(':')[0].replace('|', '/').split('/'):
                if a in ('0', '1'):
                    counts[(pop[s], li)][int(a)] += 1
        li += 1
    lp = 0.0
    for c in counts.values():
        A, N = 2, c[0] + c[1]
        lp += math.lgamma(A * ALPHA) - math.lgamma(A * ALPHA + N)
        for n in c:
            lp += math.lgamma(ALPHA + n) - math.lgamma(ALPHA)
    return lp


def main():
    ba3, vcf, meta = sys.argv[1], sys.argv[2], sys.argv[3]
    out = subprocess.run([ba3, '-V', vcf, '-M', meta, '-c', '-i', '200', '-b', '50',
                          '-n', '25', '-o', '/tmp/_collapse_check.txt'],
                         capture_output=True, text=True).stdout
    m = re.search(r'collapsed init logL.*?:\s*(-?[0-9.]+)', out)
    if not m:
        sys.exit("FAIL: BA3 did not report a collapsed init logL (is --collapse wired up?)")
    ba3_val = float(m.group(1))
    ref = independent_dirmult(vcf, meta)
    diff = abs(ba3_val - ref)
    print("BA3 collapsed init logL : %.4f" % ba3_val)
    print("independent computation : %.4f" % ref)
    print("|diff| = %.6f  ->  %s" % (diff, "PASS" if diff < 0.01 else "FAIL"))
    sys.exit(0 if diff < 0.01 else 1)


if __name__ == "__main__":
    main()
