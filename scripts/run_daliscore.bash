#!/bin/bash -e

cd ../data

for name in 1aab_1j46 BB11001
do
	reseek \
	  -daliscore_msa $name.afa \
	  -input $name.files \
	  -output ../output/daliscore.$name.tsv \
	  -log ../output/daliscore.$name.log
done
