#!/bin/bash -e

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" -output /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

reseek -msta_score 1aab_1j46.afa -input 1aab_1j46.files -output ../output/msta_score_1aab_1j46.tsv
reseek -msta_score 12.afa -input 12.files -output ../output/msta_score_12.tsv
reseek -msta_score 21.afa -input 21.files -output ../output/msta_score_21.tsv
reseek -msta_score 45.afa -input 45.files -output ../output/msta_score_45.tsv
reseek -msta_score 54.afa -input 54.files -output ../output/msta_score_54.tsv
reseek -msta_score BB11001.afa -input BB11001.files -output ../output/msta_score_BB11001.tsv
