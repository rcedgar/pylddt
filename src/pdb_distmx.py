#!/usr/bin/python3

import sys
import math
import argparse
from ms import *

Usage = \
(
"Calculate distance matrix between C-alpha atom in PDB file."
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("pdbfile", help="PDB file with structure")
Args = AP.parse_args()

scorer = create_scorer(Args)

seq, xs, ys, zs = scorer.read_pdb(Args.pdbfile)

L = len(seq)
for i in range(L):
	aa = seq[i]
	s = "%d" % i
	s += "\t" + aa
	for j in range(0, i):
		d = scorer.get_dist(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j])
		s += "\t%.1f" % d
	print(s)
