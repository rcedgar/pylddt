#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../src/lddt_model.py --ref 1.pdb --model 2.pdb \
  > ../output/lddt_model_12.txt

../src/lddt_model.py --ref 2.pdb --model 1.pdb \
  > ../output/lddt_model_21.txt

../src/lddt_model.py --ref 4.pdb --model 5.pdb \
  > ../output/lddt_model_45.txt

../src/lddt_model.py --ref 5.pdb --model 4.pdb \
  > ../output/lddt_model_54.txt

grep LDDT ../output/lddt_model*.txt
