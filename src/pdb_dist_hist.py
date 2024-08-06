#!/usr/bin/python3

import sys
import math
import argparse
from ms import *

Usage = \
(
"Calculate distance distribution for C-alpha atoms in PDB file."
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("pdbfile", help="PDB file with structure")
Args = AP.parse_args()

scorer = create_scorer(Args)

seq, xs, ys, zs = scorer.read_pdb(Args.pdbfile)

Counts = [0]*101

L = len(seq)
ds = []
for i in range(L):
	aa = seq[i]
	s = "%d" % i
	s += "\t" + aa
	for j in range(0, i):
		d = scorer.get_dist(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j])
		ds.append(d)

mind = min(ds)
maxd = max(ds)
sys.stderr.write("Distance range min %.1f to max %.1f\n" % (mind, maxd))

N = round(maxd) + 1
counts = [0]*N
for d in ds:
	counts[round(d)] += 1
print("dist\tn")
for i in range(N):
	print("%d\t%d" % (i, counts[i]))