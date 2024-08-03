#!/usr/bin/python3

import sys
import math

paths_fn = sys.argv[1]
test_msa_fn = sys.argv[2]

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

	if 0: # CA-CA bond lengths
		for i in range(0, len(xs)-1):
			x1 = xs[i-1]
			y1 = ys[i-1]
			z1 = zs[i-1]
			x = xs[i]
			y = ys[i]
			z = zs[i]
			print(get_dist(x1, y1, z1, x, y, z))
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
	D = 20
	x = 1.0/(D*D)
	w = math.exp(-x*y*y)
	return w

def DALI_dpscorefun(d1, d2):
	d0 = 0.2
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
			
pdb_fns = []
fn2data = {}
for line in open(paths_fn):
	fn = line.strip()
	pdb_fns.append(fn)
	fn2data[fn] = read_pdb(fn)
nr_pdbs = len(pdb_fns)

test_msa = read_fasta(test_msa_fn)
labels = list(test_msa.keys())
nr_cols = len(test_msa[labels[0]])
nr_test_seqs = len(labels)
msai_dx2pdb_idx = []
nr_matched = 0
msa_idxs = []
pdb_idxs = []
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
		pdb_seq, xs, ys, zs = fn2data[pdb_fn]
		if pdb_seq == test_seq:
			matched_pdb_idx = pdb_idx
			break
	if not matched_pdb_idx is None:
		msa_idxs.append(msa_idx)
		pdb_idxs.append(matched_pdb_idx)

nr_matched = len(msa_idxs)
sys.stderr.write("%d / %d MSA sequences matched to structures\n" %
				 (nr_matched, nr_test_seqs))

for i in range(nr_matched):
	msa_idxi = msa_idxs[i]
	pdb_idxi = pdb_idxs[i]
	labeli = labels[msa_idxi]
	rowi = test_msa[labeli]
	pdb_fni = pdb_fns[pdb_idxi]
	pdb_seqi, xis, yis, zis = fn2data[pdb_fni]
	for j in range(i+1, nr_matched):
		msa_idxj = msa_idxs[j]
		pdb_idxj = pdb_idxs[j]
		labelj = labels[msa_idxj]
		rowj = test_msa[labelj]
		pdb_fnj = pdb_fns[pdb_idxj]
		pdb_seqj, xjs, yjs, zjs = fn2data[pdb_fnj]

		posis, posjs = align_pair(rowi, rowj)
		dali = dali_score(posis, xis, yis, zis, posjs, xjs, yjs, zjs)
		print("%s\t%s\tdali\t%.3g" % (labeli, labelj, dali))
