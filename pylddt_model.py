#!/usr/bin/python3

import sys
import math
import argparse

Usage = \
(
'Calculate LDDT for a model structure by comparison with a reference structure'
)

AP = argparse.ArgumentParser(description = Usage)
AP.add_argument("--model", required=True, help="PDB file of predicted structure")
AP.add_argument("--ref", required=True, help="PDB file of reference structure")
AP.add_argument("--radius", required=False, type=float, default=15.0, help="Inclusion radius (R0 parameter, default 15)")
AP.add_argument("--dists", required=False, default="0.5,1,2,4", help="Distance thresholds, comma-separated (default 0.5,1,2,4)")
AP.add_argument("--output", required=False, default="/dev/stdout", help="output file (default stdout)")
Args = AP.parse_args()

model_fn = Args.model
ref_fn = Args.ref
out_fn = Args.output

R0 = Args.radius
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

def get_dist(x1, y1, z1, x2, y2, z2):
	dx = x1 - x2
	dy = y1 - y2
	dz = z1 - z2
	d2 = dx*dx + dy*dy + dz*dz
	d = math.sqrt(d2)
	return d

def lddt_score_t(pos1s, pos2s, D1, D2, t):
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
			if d1 > R0 or d2 > R0:
				continue
			nr_considered += 1
			diff = abs(d1 - d2)
			if diff < t:
				nr_preserved += 1
	fract = nr_preserved/nr_considered
	return fract

'''
Coverage: 1 (142 out of 142 residues)
Inclusion Radius: 15
Sequence separation: 0
Cutoffs: 0.5, 1, 2, 4
Chain   ResName ResNum  Asses.  Score      (Conserved/Total, over 4 thresholds)
A       GLY     0       Yes     0.617647   (42/68)
A       PHE     1       Yes     0.32       (32/100)
A       VAL     2       Yes     0.565789   (43/76)

nr_cols = 142
G     1 0.7692 40/52
F     2 0.4583 33/72
V     3 0.6053 46/76
'''

def lddt_score(seq1, pos1s, x1s, y1s, z1s, seq2, pos2s, x2s, y2s, z2s):
	nr_cols = len(pos1s)
	assert len(pos2s) == nr_cols
	print("nr_cols =", nr_cols)
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
	for coli in range(nr_cols):
		pos1i = pos1s[coli]
		pos2i = pos2s[coli]
		nr_considered = 0
		nr_preserved = 0
		for colj in range(nr_cols):
			pos1j = pos1s[colj]
			pos2j = pos2s[colj]
			for t in thresholds: # [ 0.5, 1, 2, 4 ]:
				d1 = D1[pos1i][pos1j]
				d2 = D2[pos2i][pos2j]
				if d1 > R0:
					continue
				nr_considered += 1
				diff = abs(d1 - d2)
				if diff <= t:
					nr_preserved += 1
		score = nr_preserved/nr_considered
		aa = seq1[pos1i]
		assert seq2[pos2i] == aa
		s = aa
		s += "%6d" % coli
		s += " %.4f" % score
		s += " %d/%d" % (nr_preserved, nr_considered)
		print(s)
		
		total += score
				
	avg = total/nr_cols
	return avg
			
seqi, xis, yis, zis = read_pdb(model_fn)
seqj, xjs, yjs, zjs = read_pdb(ref_fn)
if seqi != seqj:
	sys.stderr.write("\n==ERROR== PDB files have different sequences\n")
	sys.exit(1)

L = len(seqi)

posis = range(L)
posjs = range(L)

LDDT = lddt_score(seqi, posis, xis, yis, zis, seqj, posjs, xjs, yjs, zjs)
print("LDDT=%.4f" % LDDT)
