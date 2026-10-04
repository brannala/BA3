#!/usr/bin/env python3
"""Make a missing-data variant of a BA3-format input file (deterministic).

Each genotype line is made fully missing ("0 0") with probability p_full, and
half-missing ("0 a") with probability p_half, using a fixed seed so the file is
reproducible. BA3 treats a half-missing genotype as wholly missing.

usage: tests/make_missing.py in.txt out.txt [p_full=0.10] [p_half=0.02] [seed=1]
"""
import random, sys

src, dst = sys.argv[1], sys.argv[2]
p_full = float(sys.argv[3]) if len(sys.argv) > 3 else 0.10
p_half = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
with open(src) as f, open(dst, "w") as g:
    for line in f:
        t = line.split()
        if len(t) == 5:
            u = rng.random()
            if u < p_full:
                t[3] = t[4] = "0"
            elif u < p_full + p_half:
                t[3] = "0"
        g.write("   ".join(t) + "\n")
