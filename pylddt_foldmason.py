#!/usr/bin/python3

import sys
import math
import argparse

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

three2one={ 'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
			'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
			'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
			'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
			'MSE':'M' } # MSE translated to M, other special codes ignored

def read_pdb(fn):
	seq = ""
	xs = []
	ys = []
	zs = []
	for line in open(fn):
		if line.startswith("ENDMDL") or line.startswith("TER "):
			break
		name = line[12:16].strip()
		if not name == "CA":
			continue
		three = line[17:20]
		try:
			one = three2one[three]
		except:
			continue
		try:
			x = float(line[30:38])
			y = float(line[38:46])
			z = float(line[47:54])
		except:
			continue
		
		seq += one
		xs.append(x)
		ys.append(y)
		zs.append(z)

	return seq, xs, ys, zs

def read_fasta(fn):
	seqs = {}
	label = None
	seq = ""
	for line in open(fn):
		line = line.strip()
		if len(line) == 0:
			continue
		if line[0] == ">":
			if len(seq) > 0:
				seqs[label] = seq
			label = line[1:]
			seq = ""
		else:
			seq += line.replace(" ", "")
	if len(seq) > 0:
		seqs[label] = seq.upper()
	return seqs

def dali_Z(score, len1, len2):
	n12 = math.sqrt(len1*len2)
	x = min(n12, 400)
	mean = 7.9494 + 0.70852*x + 2.5895e-4*x*x - 1.9156e-6*x*x*x
	if n12 > 400:
		mean += n12 - 400
	sigma = 0.5*mean
	Z = (score - mean) / max(1.0, sigma)
	return Z

def align_pair(row1, row2):
	pos1 = 0
	pos2 = 0
	pos1s = []
	pos2s = []
	for col in range(nr_cols):
		c1 = row1[col]
		c2 = row2[col]
		gap1 = (c1 == '-' or c1 == '.')
		gap2 = (c2 == '-' or c2 == '.')
		if not gap1 and not gap2:
			pos1s.append(pos1)
			pos2s.append(pos2)
		if not gap1:
			pos1 += 1
		if not gap2:
			pos2 += 1
	return pos1s, pos2s

def get_dist(x1, y1, z1, x2, y2, z2):
	dx = x1 - x2
	dy = y1 - y2
	dz = z1 - z2
	d2 = dx*dx + dy*dy + dz*dz
	d = math.sqrt(d2)
	return d

pdb_fns = []
fn2data = {}
for line in open(paths_fn):
	fn = line.strip()
	pdb_fns.append(fn)
	fn2data[fn] = read_pdb(fn)
nr_pdbs = len(pdb_fns)

msa = read_fasta(msa_fn)
labels = list(msa.keys())
nr_cols = len(msa[labels[0]])
nr_seqs = len(labels)
msai_dx2pdb_idx = []
nr_matched = 0
msa_idxs = []
pdb_idxs = []
matched_pdb_fns = []
for msa_idx in range(nr_seqs):
	label = labels[msa_idx]
	n = len(msa[label])
	if n != nr_cols:
		assert False, "\n===ERROR=== Test MSA is not aligned\n"
	test_row = msa[label]
	test_seq = test_row.replace("-", "").replace(".", "")
	matched_pdb_idx = None
	for pdb_idx in range(nr_pdbs):
		pdb_fn = pdb_fns[pdb_idx]
		if pdb_fn in matched_pdb_fns:
			continue
		pdb_seq, xs, ys, zs = fn2data[pdb_fn]
		if pdb_seq == test_seq:
			matched_pdb_idx = pdb_idx
			matched_pdb_fns.append(pdb_fn)
			break
	if matched_pdb_idx is None:
		sys.stderr.write("\n\n===ERROR\nMSA sequence not matched to structure >%s\n\n" % label)
		sys.exit(1)
	out.write("label=%s\tpdbfile=%s\n" % (label, pdb_fn))
	msa_idxs.append(msa_idx)
	pdb_idxs.append(matched_pdb_idx)

nr_matched = len(msa_idxs)
assert nr_matched == nr_seqs

sys.stderr.write("%d / %d MSA sequences matched to structures\n" %
				 (nr_matched, nr_seqs))
if nr_matched == 0:
	sys.exit(1)

def calc_dist_mx(xs, ys, zs):
	n = len(xs)
	mx = []
	for i in range(n):
		mx.append([])
		for j in range(n):
			mx[i].append(0)

	assert len(ys) == n
	assert len(zs) == n
	for i in range(n):
		for j in range(i+1, n):
			d = get_dist(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j]);
			mx[i][j] = d
			mx[j][i] = d
	return mx

def calc_col2pos(row):
	col2pos = []
	pos = 0
	for c in row:
		if c == '-' or c == '.':
			col2pos.append(None)
		else:
			col2pos.append(pos)
			pos += 1
	return col2pos

dist_mxs= []
col2pos_vec = []
for i in range(nr_seqs):
	msa_idx = msa_idxs[i]
	pdb_idx = pdb_idxs[i]
	pdb_fn = pdb_fns[pdb_idx]
	pdb_seq, xs, ys, zs = fn2data[pdb_fn]
	dist_mxs.append(calc_dist_mx(xs, ys, zs))
	label = labels[msa_idx]
	row = msa[label]
	col2pos_vec.append(calc_col2pos(row))

def get_row(i):
	label = labels[i]
	return msa[label]

def calc_col_score(col):
	total_pair_scores = 0
	nr_seq_pairs = 0
	for seq_idxi in range(nr_seqs):
		posi = col2pos_vec[seq_idxi][col]
		if posi is None:
			continue
		dist_mxi = dist_mxs[seq_idxi]
		for seq_idxj in range(seq_idxi+1, nr_seqs):
			posj = col2pos_vec[seq_idxj][col]
			if posj is None:
				continue
			nr_seq_pairs += 1
			dist_mxj = dist_mxs[seq_idxj]
			pair_score = 0

			nr_pairs = 0
			score = 0
			for col2 in range(nr_cols):
				if col2 == col:
					continue
				posi2 = col2pos_vec[seq_idxi][col2]
				posj2 = col2pos_vec[seq_idxj][col2]
				if posi2 is None or posj2 is None:
					continue
				
				di = dist_mxi[posi][posi2]
				dj = dist_mxj[posj][posj2]
				if symmetry == "first":
					if di > R0:
						continue
				elif symmetry == "both":
					if di > R0 and dj > R0:
						continue
				elif symmetry == "either":
					if di > R0 or dj > R0:
						continue

				nr_pairs += 1
				d_l = abs(di - dj)
				score += (int(d_l < 0.5) + int(d_l < 1.0) + int(d_l < 2.0) + int(d_l < 4.0))/4

			if nr_pairs > 0:
				pair_score = score/nr_pairs
				total_pair_scores += pair_score
	if nr_seq_pairs == 0:
		return 0
	return total_pair_scores/nr_seq_pairs

total_col_scores = 0
for col in range(nr_cols):
	col_score = calc_col_score(col)
	total_col_scores += col_score
	# print(col, "%.4f" % col_score)

mean_score = total_col_scores/nr_cols
out.write("LDDT=%.4f\n" % mean_score)
out.close()
