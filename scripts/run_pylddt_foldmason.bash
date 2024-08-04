#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../pylddt_foldmason.py --msa 1aab_1j46.afa --pdbfiles 1aab_1j46.files > ../output/pylddt_foldmason_1aab_1j46.tsv
../pylddt_foldmason.py --msa 12.afa --pdbfiles 12.files > ../output/pylddt_foldmason_12.tsv
../pylddt_foldmason.py --msa 21.afa --pdbfiles 21.files > ../output/pylddt_foldmason_21.tsv
../pylddt_foldmason.py --msa 45.afa --pdbfiles 45.files > ../output/pylddt_foldmason_45.tsv
../pylddt_foldmason.py --msa 54.afa --pdbfiles 54.files > ../output/pylddt_foldmason_54.tsv
../pylddt_foldmason.py --msa BB11001.afa --pdbfiles BB11001.files > ../output/pylddt_foldmason_BB11001.tsv

grep LDDT= ../output/pylddt_foldmason*
