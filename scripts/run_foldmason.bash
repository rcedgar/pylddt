#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../scripts/foldmason_msa2lddt.bash 1aab_1j46.afa 1aab_1j46.files \
  > ../output/foldmason_1aab_1j46.txt

../scripts/foldmason_msa2lddt.bash 12.afa 12.files \
  > ../output/foldmason_12.txt

../scripts/foldmason_msa2lddt.bash 21.afa 21.files \
  > ../output/foldmason_21.txt

../scripts/foldmason_msa2lddt.bash 45.afa 45.files \
  > ../output/foldmason_45.txt

../scripts/foldmason_msa2lddt.bash 54.afa 54.files \
  > ../output/foldmason_54.txt

## ../scripts/foldmason_msa2lddt.bash BB11001.afa BB11001.files \
## > ../output/foldmason_BB11001.txt

cd ../output
grep "Average MSA LDDT" foldmason_*.txt
