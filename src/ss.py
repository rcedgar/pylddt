#!/usr/bin/python3

import sys
import math

class StructureScorer:
	def __init__(self):
		self.three2one = { 'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
					'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
					'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
					'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
					'MSE':'M' } # MSE translated to M, other special codes ignored
		self.symmetry = "first"
		self.R0 = 15
		self.dists = [ 0.5, 1, 2, 4 ]

	def read_pdb(self, fn):
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
				one = self.three2one[three]
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

	def read_msa(self, fn):
		sys.stderr.write("\n");
		self.msa = {}
		self.labels = []
		label = None
		seq = ""
		for line in open(fn):
			line = line.strip()
			if len(line) == 0:
				continue
			if line[0] == ">":
				if len(seq) > 0:
					self.labels.append(label)
					self.msa[label] = seq.upper()
				label = line[1:]
				seq = ""
			else:
				seq += line.replace(" ", "")
		if len(seq) > 0:
			self.labels.append(label)
			self.msa[label] = seq.upper()
		self.nr_seqs = len(self.labels)
		if self.nr_seqs == 0:
			sys.stderr.write("\n=== ERROR=== empty MSA\n\n")
			sys.exit(1)
		self.nr_cols = len(self.msa[self.labels[0]])
		for i in range(1, self.nr_seqs):
			if len(self.msa[self.labels[i]]) != self.nr_cols:
				sys.stderr.write("\n=== ERROR=== FASTA is not aligned\n\n")
				sys.exit(1)

	def dali_Z(score, len1, len2):
		n12 = math.sqrt(len1*len2)
		x = min(n12, 400)
		mean = 7.9494 + 0.70852*x + 2.5895e-4*x*x - 1.9156e-6*x*x*x
		if n12 > 400:
			mean += n12 - 400
		sigma = 0.5*mean
		Z = (score - mean) / max(1.0, sigma)
		return Z

	def align_pair(self, row1, row2):
		pos1 = 0
		pos2 = 0
		pos1s = []
		pos2s = []
		for col in range(self, self.nr_cols):
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

	def get_dist(self, x1, y1, z1, x2, y2, z2):
		dx = x1 - x2
		dy = y1 - y2
		dz = z1 - z2
		d2 = dx*dx + dy*dy + dz*dz
		d = math.sqrt(d2)
		return d

	def calc_dist_mx(self, xs, ys, zs):
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
				d = self.get_dist(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j]);
				mx[i][j] = d
				mx[j][i] = d
		return mx

	def calc_col2pos(self, row):
		col2pos = []
		pos = 0
		for c in row:
			if c == '-' or c == '.':
				col2pos.append(None)
			else:
				col2pos.append(pos)
				pos += 1
		return col2pos

	def read_pdbs(self, paths_fn):
		self.pdb_fns = []
		self.fn2data = {}
		for line in open(paths_fn):
			fn = line.strip()
			self.pdb_fns.append(fn)
			self.fn2data[fn] = self.read_pdb(fn)
		self.nr_pdbs = len(self.pdb_fns)
	
	def match_seqs_to_pdbs(self):
		self.msai_dx2pdb_idx = []
		self.nr_matched = 0
		self.msa_idxs = []
		self.pdb_idxs = []
		self.matched_pdb_fns = []
		self.not_matched_labels = []
		for msa_idx in range(self.nr_seqs):
			label = self.labels[msa_idx]
			n = len(self.msa[label])
			if n != self.nr_cols:
				sys.stderr.write("\n===ERROR=== Test MSA is not aligned\n")
				sys.exit(1)
			test_row = self.msa[label]
			test_seq = test_row.replace("-", "").replace(".", "")
			matched_pdb_idx = None
			for pdb_idx in range(self.nr_pdbs):
				pdb_fn = self.pdb_fns[pdb_idx]
				if pdb_fn in self.matched_pdb_fns:
					continue
				pdb_seq, _, _, _s = self.fn2data[pdb_fn]
				if pdb_seq == test_seq:
					matched_pdb_idx = pdb_idx
					self.matched_pdb_fns.append(pdb_fn)
					self.nr_matched += 1
					break
			if matched_pdb_idx is None:
				self.not_matched_labels.append(label)
				continue
			self.msa_idxs.append(msa_idx)
			self.pdb_idxs.append(matched_pdb_idx)
		if self.nr_matched < 2:
			sys.stderr.write("\n===ERROR=== %d/%d sequences matched, must be >=2\n" %
					(self.nr_matched, self.nr_seqs))
			sys.exit(1)

	def set_dist_mxs(self):
		self.dist_mxs = []
		for i in range(self.nr_seqs):
			pdb_idx = self.pdb_idxs[i]
			pdb_fn = self.pdb_fns[pdb_idx]
			pdb_seq, xs, ys, zs = self.fn2data[pdb_fn]
			self.dist_mxs.append(self.calc_dist_mx(xs, ys, zs))

	def set_col2pos_vec(self):
		self.col2pos_vec = []
		for i in range(self.nr_seqs):
			msa_idx = self.msa_idxs[i]
			label = self.labels[msa_idx]
			row = self.msa[label]
			self.col2pos_vec.append(self.calc_col2pos(row))

	def calc_col_score(self, col):
		total_pair_scores = 0
		nr_seq_pairs = 0
		for seq_idxi in range(self.nr_seqs):
			posi = self.col2pos_vec[seq_idxi][col]
			if posi is None:
				continue
			dist_mxi = self.dist_mxs[seq_idxi]
			for seq_idxj in range(seq_idxi+1, self.nr_seqs):
				posj = self.col2pos_vec[seq_idxj][col]
				if posj is None:
					continue
				nr_seq_pairs += 1
				dist_mxj = self.dist_mxs[seq_idxj]
				pair_score = 0

				nr_pairs = 0
				score = 0
				for col2 in range(self.nr_cols):
					if col2 == col:
						continue
					posi2 = self.col2pos_vec[seq_idxi][col2]
					posj2 = self.col2pos_vec[seq_idxj][col2]
					if posi2 is None or posj2 is None:
						continue
				
					di = dist_mxi[posi][posi2]
					dj = dist_mxj[posj][posj2]
					if self.symmetry == "first":
						if di > self.R0:
							continue
					elif self.symmetry == "both":
						if di > self.R0 and dj > self.R0:
							continue
					elif self.symmetry == "either":
						if di > self.R0 or dj > self.R0:
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
	
	def calc_mean_col_score(self):
		self.total_col_scores = 0
		self.col_scores = []
		for col in range(self.nr_cols):
			col_score = self.calc_col_score(col)
			self.col_scores.append(col_score)
			self.total_col_scores += col_score
		self.mean_col_score = self.total_col_scores/self.nr_cols
		return self.mean_col_score
