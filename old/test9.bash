#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../src/daliz.py --msa 1aab_1j46.afa --pdbfiles 1aab_1j46.files
