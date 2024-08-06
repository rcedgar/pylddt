#!/usr/bin/python3

import sys
import math
import argparse
from ms import *

Usage = \
(
'Calculate LDDT for a sequence alignment by comparison with '
'reference structures as average of pair-wise LDDT.'
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--msa", required=True, help="Multiple Sequence Alignment (FASTA format)")
AP.add_argument("--pdbfiles", required=True, help="Text file with one PDB pathname per line")
AP.add_argument("--radius", required=False, type=float, default=15.0, help="LDDT inclusion radius (default 15)")
AP.add_argument("--dists", required=False, default="0.5,1,2,4", help="LDDT distance thresholds, comma-separated (default 0.5,1,2,4)")
AP.add_argument("--symmetry", required=False, choices=[ "first", "both", "either" ], default="first", help="Set R0 according to both / either / first (default first)")
AP.add_argument("--pairs", action="store_true", help="Report pair-wise LDDTs (default don't show)")
Args = AP.parse_args()

msa_fn = Args.msa
paths_fn = Args.pdbfiles

scorer = create_scorer(Args)

scorer.read_pdbs(paths_fn)
scorer.read_msa(msa_fn)
scorer.match_seqs_to_pdbs()

nr_matched = len(scorer.msa_idxs)
if nr_matched != scorer.nr_seqs:
	sys.stderr.write("\n===ERROR=== %d / %d MSA sequences matched to structures\n" %
					 (scorer.nr_matched, scorer.nr_seqs))
	sys.exit(1)

scorer.set_dist_mxs()
scorer.set_col2pos_vec()
scorer.match_seqs_to_pdbs()
scorer.calc_lddt_scores()

if Args.pairs:
	n = len(scorer.LDDT_label1s)
	for i in range(n):
		print("%s\t%s\t%.4f" %
			 (scorer.LDDT_label1s[i], scorer.LDDT_label2s[i], scorer.LDDT_scores[i]))

print("LDDT_rce_mean\t%.4f\t%s\n" % (scorer.mean_LDDT_score, msa_fn))
