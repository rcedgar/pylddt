#!/usr/bin/python3

import sys

qpos1s = []
qpos2s = []
tpos1s = []
tpos2s = []
distances = []
dist_subs = []
d_ls = []
scores = []
includes = []

nr_cols = None
for line in open(sys.argv[1]):
	line = line.strip()
	line = line.replace(",", "")
	if line.startswith("alnLength="):
		nr_cols = int(line.replace("alnLength=", ""))
		continue
	if not line.startswith("query_idx1="):
		continue
	flds = line.split(" ")
	qpos1 = None
	qpos2 = None
	tpos1 = None
	tpos2 = None
	distance = None
	dist_sub = None
	d_l = None
	score = None
	include = None
	for fld in flds:
		if fld.startswith("query_idx1="):
			qpos1 = int(fld.replace("query_idx1=", ""))
		if fld.startswith("query_idx2="):
			qpos2 = int(fld.replace("query_idx2=", ""))
		if fld.startswith("tpos1="):
			tpos1 = int(fld.replace("tpos1=", ""))
		if fld.startswith("tpos2="):
			tpos2 = int(fld.replace("tpos2=", ""))
		if fld.startswith("distance="):
			distance = float(fld.replace("distance=", ""))
		if fld.startswith("dist_sub="):
			dist_sub = float(fld.replace("dist_sub=", ""))
		if fld.startswith("score="):
			score = float(fld.replace("score=", ""))
		if fld.startswith("d_l="):
			d_l= float(fld.replace("d_l=", ""))
		if fld.startswith("include="):
			include = fld.replace("include=", "")

	qpos1s.append(qpos1)
	qpos2s.append(qpos2)
	tpos1s.append(tpos1)
	tpos2s.append(tpos2)
	distances.append(distance)
	dist_subs.append(dist_sub)
	d_ls.append(d_l)
	scores.append(score)
	includes.append(include)

N = len(qpos1s)

include_mx = []
tpos1_mx = []
tpos2_mx = []
distance_mx = []
dist_sub_mx = []
d_l_mx = []
score_mx = []
for col in range(nr_cols):
	include_mx.append([])
	tpos1_mx.append([])
	tpos2_mx.append([])
	distance_mx.append([])
	dist_sub_mx.append([])
	d_l_mx.append([])
	score_mx.append([])
	for col2 in range(nr_cols):
		include_mx[col].append(None)
		tpos1_mx[col].append(None)
		tpos2_mx[col].append(None)
		distance_mx[col].append(None)
		dist_sub_mx[col].append(None)
		d_l_mx[col].append(None)
		score_mx[col].append(None)

for i in range(N):
	qpos1 = qpos1s[i]
	qpos2 = qpos2s[i]
	include = includes[i]
	tpos1 = tpos1s[i]
	tpos2 = tpos2s[i]
	include = includes[i]
	dist_sub = dist_subs[i]
	d_l = d_ls[i]
	score = scores[i]

	include_mx[qpos1][qpos2] = include
	tpos1_mx[qpos1][qpos2] = tpos1
	tpos2_mx[qpos1][qpos2] = tpos2
	distance_mx[qpos1][qpos2] = distance
	dist_sub_mx[qpos1][qpos2] = dist_sub
	d_l_mx[qpos1][qpos2] = d_l
	score_mx[qpos1][qpos2] = score

	include_mx[qpos2][qpos1] = include
	tpos1_mx[qpos2][qpos1] = tpos1
	tpos2_mx[qpos2][qpos1] = tpos2
	distance_mx[qpos2][qpos1] = distance
	dist_sub_mx[qpos2][qpos1] = dist_sub
	d_l_mx[qpos2][qpos1] = d_l
	score_mx[qpos2][qpos1] = score

for qpos1 in range(nr_cols):
	for qpos2 in range(qpos1+1, nr_cols):
		include = include_mx[qpos1][qpos2]
		tpos1 = tpos1_mx[qpos1][qpos2]
		tpos2 = tpos2_mx[qpos1][qpos2]
		distance = distance_mx[qpos1][qpos2]
		dist_sub = dist_sub_mx[qpos1][qpos2]
		d_l = d_l_mx[qpos1][qpos2]
		score = score_mx[qpos1][qpos2]

		s = "qpos1=%d" % qpos1
		s += " qpos2=%d" % qpos2
		if include is None:
			s += " include=None"
		else:
			s += " include=%c" % include
		if include == "T":
			s += " tpos1=%s" % str(tpos1)
			s += " tpos2=%s" % str(tpos2)
			s += " distance=%s" % str(distance)
			s += " dist_sub=%s" % str(dist_sub)
			s += " d_l=%s" % str(d_l)
			s += " score=%s" % str(score)
		print(s)
