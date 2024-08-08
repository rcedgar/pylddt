#!/usr/bin/python3

import sys
import math
import argparse
from ms import *

Usage = \
(
"Calculate DALI Z score from a test Multiple Sequence Alignment (MSA) and PDB files with structures. "
"Reports Z for each pair of sequences, and the mean Z over all pairs. "
"Sequences are matched to structures by string matching, the sequence of a PDB file is taken from CA ATOM records. "
"Structures which do not match the MSA are ignored. "
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--msa", required=True, help="Multiple Sequence Alignment (FASTA format)")
AP.add_argument("--pdbfiles", required=True, help="Text file with one PDB pathname per line")
AP.add_argument("--radius", required=False, type=float, default=None, help="Inclusion radius as for LDDT (default none)")
AP.add_argument("--symmetry", required=False, choices=[ "first", "both", "either" ], default="first", 
                help="Which distance(s) considered for radius (default first)")
AP.add_argument("--horizon", required=False, type=float, default=20.0, help="DALI horizon (default 20)")
AP.add_argument("--diagwt", required=False, type=float, default=0.2, help="DALI diagonal weight (default 0.2)")
Args = AP.parse_args()

scorer = create_scorer(Args)

scorer.read_pdbs(Args.pdbfiles)
scorer.read_msa(Args.msa)
scorer.match_seqs_to_pdbs()
for label in scorer.not_matched_labels:
    print("MSA_label_not_matched=%s" % label)

scorer.calc_dali_scores()
n = len(scorer.DALI_label1s)
for i in range(n):
    s = "label1=" + scorer.DALI_label1s[i]
    s += "\tlabel2=" + scorer.DALI_label2s[i]
    s += "\tscore=%.0f" % scorer.DALI_scores[i]
    s += "\tZ=%.3f" % scorer.DALI_Zs[i]
    print(s)

print("MSA=%s\tmean_DALI_score=%.0f\tmean_DALI_Z=%.3f" %
      (Args.msa, scorer.mean_DALI_score, scorer.mean_DALI_Z))
