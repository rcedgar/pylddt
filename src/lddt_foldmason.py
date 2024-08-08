#!/usr/bin/python3

import sys
import math
import argparse
from ms import *

Usage = \
(
'Calculate LDDT for a sequence alignment by comparison with '
'reference structures, using the Foldmason method.'
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--msa", required=True, help="Multiple Sequence Alignment (FASTA format)")
AP.add_argument("--pdbfiles", required=True, help="Text file with one PDB pathname per line")
AP.add_argument("--radius", required=False, type=float, default=15.0, help="LDDT inclusion radius (default 15)")
AP.add_argument("--dists", required=False, default="0.5,1,2,4", help="LDDT distance thresholds, comma-separated (default 0.5,1,2,4)")
AP.add_argument("--symmetry", required=False, choices=[ "first", "both", "either" ], default="first", help="Set R0 according to both / either / first (default first)")
AP.add_argument("--cols", action="store_true", help="Report column scores")
Args = AP.parse_args()

msa_fn = Args.msa
paths_fn = Args.pdbfiles

scorer = create_scorer(Args)

scorer.read_pdbs(paths_fn)
scorer.read_msa(msa_fn)
scorer.match_seqs_to_pdbs()

nr_matched = len(scorer.msa_idxs)
if nr_matched != scorer.nr_seqs:
	sys.stderr.write("\n===ERROR===\n%d / %d MSA sequences matched to structures\n" %
					 (scorer.nr_matched, scorer.nr_seqs))
	sys.stderr.write("first not matched >%s\n" %
					 (scorer.not_matched_labels[0]))
	sys.exit(1)

scorer.set_dist_mxs()
scorer.set_col2pos_vec()
LDDT = scorer.calc_mean_col_score()

if Args.cols:
	print("col_nr\tcol_str\tscore")
	for col_nr in range(scorer.nr_cols):
		col_str = scorer.msa_col(col_nr)
		print("%d\t%s\t%.4f" % (col_nr, col_str, scorer.col_scores[col_nr]))

print("LDDT_foldmason\t%.4f\t%s\n" % (LDDT, msa_fn))
