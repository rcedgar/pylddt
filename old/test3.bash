#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../src/lddt_foldmason.py --msa BB11001.afa --pdbfiles BB11001.files --cols yes
