#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../src/lddt_foldmason.py --msa 1aab_1j46.afa --pdbfiles 1aab_1j46.files > ../output/lddt_foldmason_1aab_1j46.tsv
../src/lddt_foldmason.py --msa 12.afa --pdbfiles 12.files > ../output/lddt_foldmason_12.tsv
../src/lddt_foldmason.py --msa 21.afa --pdbfiles 21.files > ../output/lddt_foldmason_21.tsv
../src/lddt_foldmason.py --msa 45.afa --pdbfiles 45.files > ../output/lddt_foldmason_45.tsv
../src/lddt_foldmason.py --msa 54.afa --pdbfiles 54.files > ../output/lddt_foldmason_54.tsv
../src/lddt_foldmason.py --msa BB11001.afa --pdbfiles BB11001.files > ../output/lddt_foldmason_BB11001.tsv

grep LDDT ../output/lddt_foldmason*.tsv
