#!/usr/bin/python3

import sys
import math
import argparse
from lddt import LDDT

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
AP.add_argument("--output", required=False, default="/dev/stdout", help="output file (default stdout)")
Args = AP.parse_args()

msa_fn = Args.msa
paths_fn = Args.pdbfiles
out_fn = Args.output
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

out = open(out_fn, "w")

ld = LDDT()
ld.read_pdbs(paths_fn)
ld.read_msa(msa_fn)
ld.match_seqs_to_pdbs()

nr_matched = len(ld.msa_idxs)
assert nr_matched == ld.nr_seqs

sys.stderr.write("%d / %d MSA sequences matched to structures\n" %
				 (ld.nr_matched, ld.nr_seqs))
if nr_matched == 0:
	sys.exit(1)

ld.set_dist_mxs()
ld.set_col2pos_vec()
LDDT = ld.calc_mean_col_score()

out.write("LDDT=%.4f\n" % LDDT)
out.close()
