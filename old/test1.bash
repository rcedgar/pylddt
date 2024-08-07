#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../src/lddt_foldmason.py --msa 12.afa --pdbfiles 12.files
../src/lddt_rce.py --msa 12.afa --pdbfiles 12.files
