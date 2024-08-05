#!/usr/bin/python3

import sys
import math
import argparse
from ms import *

Usage = \
(
'Calculate LDDT for a model structure by comparison with a reference structure'
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--ref", required=True, help="PDB file of reference structure")
AP.add_argument("--model", required=True, help="PDB file of predicted structure")
AP.add_argument("--radius", required=False, type=float, default=15.0, help="Inclusion radius (R0 parameter, default 15)")
AP.add_argument("--dists", required=False, default="0.5,1,2,4", help="Distance thresholds, comma-separated (default 0.5,1,2,4)")
AP.add_argument("--output", required=False, default="/dev/stdout", help="output file (default stdout)")
AP.add_argument("--cols", required=False, default=False, choices=[ "yes", "no" ],  help="Report column scores yes/no (default no)")

Args = AP.parse_args()

model_fn = Args.model
ref_fn = Args.ref
out_fn = Args.output

# For standard LDDT, R0 is defined by the reference
Args.symmetry = "first"
scorer = create_scorer(Args)

seqi, xis, yis, zis = scorer.read_pdb(ref_fn)
seqj, xjs, yjs, zjs = scorer.read_pdb(model_fn)

L = len(seqi)
if len(seqj) != L:
    sys.stderr.write("\n===ERROR=== sequences different lengths\n")
    sys.exit(1)

posis = range(L)
posjs = range(L)

LDDT = scorer.lddt_score(seqi, posis, xis, yis, zis, seqj, posjs, xjs, yjs, zjs)
if Args.cols == "yes":
	print("col\taa\tnr_pres\tnr_cons\tscore")
	for i in range(scorer.nr_cols):
		c = seqi[i]
		print("%d\t%c\t%d\t%d\t%.4f" % (i, c, scorer.nr_preserveds[i], scorer.nr_considereds[i], scorer.col_scores[i]))

print("LDDT\t%.4f" % LDDT)
