#!/usr/bin/python3

import sys
import math
import argparse
from ss import StructureScorer

Usage = \
(
"Usage TODO."
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--msa", required=True, help="Multiple Sequence Alignment (FASTA format)")
AP.add_argument("--pdbfiles", required=True, help="Text file with one PDB pathname per line")
AP.add_argument("--radius", required=False, type=float, default=15.0, help="LDDT inclusion radius (default 15)")
AP.add_argument("--dists", required=False, default="0.5,1,2,4", help="LDDT distance thresholds, comma-separated (default 0.5,1,2,4)")
AP.add_argument("--symmetry", required=False, choices=[ "first", "both", "either" ], default="first", help="Set R0/dali_radius according to both / either / first (default first)")
AP.add_argument("--cols", required=False, default=False, choices=[ "yes", "no" ],  help="Report column scores yes/no (default no)")
Args = AP.parse_args()

msa_fn = Args.msa
paths_fn = Args.pdbfiles
R0 = Args.radius
dists = Args.dists
symmetry = Args.symmetry
thresholds = []
flds = dists.split(',')
for fld in flds:
	try:
		t = float(fld)
	except:
		sys.stderr.write("Invalid distance in --dists argument\n")
		sys.exit(0)
	thresholds.append(t)

scorer = StructureScorer()

scorer.R0 = Args.radius
scorer.dists = Args.dists
scorer.symmetry = Args.symmetry

scorer.read_pdbs(paths_fn)
scorer.read_msa(msa_fn)
scorer.match_seqs_to_pdbs()

nr_matched = len(scorer.msa_idxs)
assert nr_matched == scorer.nr_seqs

sys.stderr.write("%d / %d MSA sequences matched to structures\n" %
				 (scorer.nr_matched, scorer.nr_seqs))
if nr_matched == 0:
	sys.exit(1)

scorer.set_dist_mxs()
scorer.set_col2pos_vec()
score = scorer.calc_mean_col_score()

if Args.cols == "yes":
	for i in range(scorer.nr_cols):
		print("%d\t%.4f" % (i, scorer.col_scores[i]))

print("LDDT=%.4f\n" % score)
