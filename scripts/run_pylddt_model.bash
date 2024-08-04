#!/bin/bash

if [ ! -d ../data ] ; then
	echo "Must run $0 from scripts/ directory" > /dev/stderr
	exit 1
fi

cd ../data
mkdir -p ../output

../pylddt_model.py --model 1.pdb --ref 2.pdb \
  > ../output/pylddt_model_12.txt

../pylddt_model.py --model 2.pdb --ref 1.pdb \
  > ../output/pylddt_model_21.txt

../pylddt_model.py --model 4.pdb --ref 5.pdb \
  > ../output/pylddt_model_45.txt

../pylddt_model.py --model 5.pdb --ref 4.pdb \
  > ../output/pylddt_model_54.txt
