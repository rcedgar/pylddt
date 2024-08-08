#!/bin/bash -e

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../src/daliz.py --msa 1aab_1j46.afa --pdbfiles 1aab_1j46.files > ../output/daliz_1aab_1j46.tsv
../src/daliz.py --msa 12.afa --pdbfiles 12.files > ../output/daliz_12.tsv
../src/daliz.py --msa 21.afa --pdbfiles 21.files > ../output/daliz_21.tsv
../src/daliz.py --msa 45.afa --pdbfiles 45.files > ../output/daliz_45.tsv
../src/daliz.py --msa 54.afa --pdbfiles 54.files > ../output/daliz_54.tsv
../src/daliz.py --msa BB11001.afa --pdbfiles BB11001.files > ../output/daliz_BB11001.tsv

grep mean_DALI_Z ../output/daliz*.tsv
