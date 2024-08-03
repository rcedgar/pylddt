#!/usr/bin/python3

import sys
import math
import argparse

Usage = \
(
"Calculate DALI Z score and LDDT from a test Multiple Sequence Alignment (MSA) and PDB files with structures. "
"Reports Z and LDDT for each pair of sequences, and the mean Z and LDDT over all pairs. "
"Sequences are matched to structures by string matching, the sequence of a PDB file is taken from CA ATOM records. "
"Structures which do not match the MSA are ignored. "
"Sequences in the MSA which do not match a structure are reported, then ignored."
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--msa", required=True, help="Multiple Sequence Alignment (FASTA format)")
AP.add_argument("--pdbfiles", required=True, help="Text file with one PDB pathname per line")
AP.add_argument("--radius", required=False, type=float, default=15.0, help="LDDT inclusion radius (default 15)")
AP.add_argument("--dists", required=False, default="0.5,1,2,4", help="LDDT distance thresholds, comma-separated (default 0.5,1,2,4)")
AP.add_argument("--horizon", required=False, type=float, default=20.0, help="DALI horizon (default 20)")
AP.add_argument("--diagwt", required=False, type=float, default=0.2, help="DALI diagonal weight (default 0.2)")
AP.add_argument("--output", required=False, default="/dev/stdout", help="output file (default stdout)")
Args = AP.parse_args()

msa_fn = Args.msa
paths_fn = Args.pdbfiles
out_fn = Args.output
R0 = Args.radius
D = Args.horizon
d0 = Args.diagwt
dists = Args.dists
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

def DALI_weight(y):
	x = 1.0/(D*D)
	w = math.exp(-x*y*y)
	return w

def DALI_dpscorefun(d1, d2):
	diff = abs(d1 - d2)
	mean = (d1 + d2)/2
	if mean > 100:
		return 0
	w = DALI_weight(mean)
	ratio = diff/mean
	if mean == 0:
		return w*d0
	return w*(d0 - ratio)

def dali_score(pos1s, x1s, y1s, z1s, pos2s, x2s, y2s, z2s):
	n = len(pos1s)
	assert len(pos2s) == n
	score = 0
	for coli in range(n):
		pos1i = pos1s[coli]
		pos2i = pos2s[coli]
		for colj in range(0, n):
			if colj == coli:
				continue
			pos1j = pos1s[colj]
			pos2j = pos2s[colj]

			dij1 = get_dist(
				x1s[pos1i], y1s[pos1i], z1s[pos1i],
				x1s[pos1j], y1s[pos1j], z1s[pos1j])

			dij2 = get_dist(
				x2s[pos2i], y2s[pos2i], z2s[pos2i],
				x2s[pos2j], y2s[pos2j], z2s[pos2j])
			
			score += DALI_dpscorefun(dij1, dij2)

	score += n*0.2
	return score

def lddt_score_threadhold(pos1s, pos2s, D1, D2, t):
	nr_cols = len(pos1s)
	assert len(pos2s) == nr_cols
	nr_considered = 0
	nr_preserved = 0
	for i in range(nr_cols):
		pos1i = pos1s[i]
		pos2i = pos2s[i]
		for j in range(i+1, nr_cols):
			pos1j = pos1s[j]
			pos2j = pos2s[j]
			d1 = D1[pos1i][pos1j]
			d2 = D2[pos2i][pos2j]
			if d1 < R0 or d2 < R0:
				continue
			nr_considered += 1
			diff = abs(d1 - d2)
			if diff <= t:
				nr_preserved += 1
	fract = nr_preserved/nr_considered
	return fract

def lddt_score(pos1s, x1s, y1s, z1s, pos2s, x2s, y2s, z2s):
	L1 = len(x1s)
	L2 = len(x2s)

	assert len(y1s) == L1
	assert len(z1s) == L1

	assert len(y2s) == L2
	assert len(z2s) == L2
	
	D1 = []
	for i in range(L1):
		D1.append([ None ]*L1)

	D2 = []
	for i in range(L2):
		D2.append([ None ]*L2)
		
	for i in range(L1):
		for j in range(i, L1):
			d = get_dist(x1s[i], y1s[i], z1s[i], x1s[j], y1s[j], z1s[j])
			D1[i][j] = d
			D1[j][i] = d

	for i in range(L2):
		for j in range(i, L2):
			d = get_dist(x2s[i], y2s[i], z2s[i], x2s[j], y2s[j], z2s[j])
			D2[i][j] = d
			D2[j][i] = d

	total = 0
	for t in thresholds: # [ 0.5, 1, 2, 4 ]:
		score = lddt_score_threadhold(pos1s, pos2s, D1, D2, t)
		total += score
	avg = total/4
	return avg
			
pdb_fns = []
fn2data = {}
for line in open(paths_fn):
	fn = line.strip()
	pdb_fns.append(fn)
	fn2data[fn] = read_pdb(fn)
nr_pdbs = len(pdb_fns)

test_msa = read_fasta(msa_fn)
labels = list(test_msa.keys())
nr_cols = len(test_msa[labels[0]])
nr_test_seqs = len(labels)
msai_dx2pdb_idx = []
nr_matched = 0
msa_idxs = []
pdb_idxs = []
matched_pdb_fns = []
for msa_idx in range(nr_test_seqs):
	label = labels[msa_idx]
	n = len(test_msa[label])
	if n != nr_cols:
		assert False, "\n===ERROR=== Test MSA is not aligned\n"
	test_row = test_msa[label]
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
		out.write("label=%s\tpdbfile=SEQUENCE_NOT_MATCHED\n" % label)
	else:
		out.write("label=%s\tpdbfile=%s\n" % (label, pdb_fn))
		msa_idxs.append(msa_idx)
		pdb_idxs.append(matched_pdb_idx)

nr_matched = len(msa_idxs)
sys.stderr.write("%d / %d MSA sequences matched to structures\n" %
				 (nr_matched, nr_test_seqs))
if nr_matched == 0:
	sys.exit(1)

total_Z = 0
total_LDDT = 0
nr_pairs = 0
for i in range(nr_matched):
	msa_idxi = msa_idxs[i]
	pdb_idxi = pdb_idxs[i]
	labeli = labels[msa_idxi]
	rowi = test_msa[labeli]
	pdb_fni = pdb_fns[pdb_idxi]
	pdb_seqi, xis, yis, zis = fn2data[pdb_fni]
	Li = len(pdb_seqi)
	for j in range(i+1, nr_matched):
		msa_idxj = msa_idxs[j]
		pdb_idxj = pdb_idxs[j]
		labelj = labels[msa_idxj]
		rowj = test_msa[labelj]
		pdb_fnj = pdb_fns[pdb_idxj]
		pdb_seqj, xjs, yjs, zjs = fn2data[pdb_fnj]
		Lj = len(pdb_seqj)
		posis, posjs = align_pair(rowi, rowj)
		score = dali_score(posis, xis, yis, zis, posjs, xjs, yjs, zjs)
		Z = dali_Z(score, Li, Lj)
		LDDT = lddt_score(posis, xis, yis, zis, posjs, xjs, yjs, zjs)
		total_Z += Z
		total_LDDT += LDDT
		nr_pairs += 1
		out.write("label1=%s\tlabel2=%s\tDALI_Z=%.1f\tLDDT=%.4f\n" %
			(labeli, labelj, Z, LDDT))

mean_Z = total_Z/nr_pairs
mean_LDDT = total_LDDT/nr_pairs

out.write("MSA=%s\tmatched=%d/%d\tpairs=%d\tmean_Z=%.1f\tmean_LDDT=%.4f" %
	  (msa_fn, nr_matched, nr_test_seqs, nr_pairs, mean_Z, mean_LDDT))
