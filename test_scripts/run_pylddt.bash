#!/bin/bash -e

if [ ! -d ../test_data ] ; then
	echo "Must run $0 from test_scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../test_data
mkdir -p ../test_output

../pylddt.py --msa 1aab_1j46.afa --pdbfiles 1aab_1j46.files > ../test_output//pylddt_1aab_1j46.tsv
../pylddt.py --msa 12.afa --pdbfiles 12.files > ../test_output/pylddt_12.tsv
../pylddt.py --msa 45.afa --pdbfiles 45.files > ../test_output/pylddt_45.tsv
../pylddt.py --msa BB11001.afa --pdbfiles BB11001.files > ../test_output/pylddt_BB11001.tsv
