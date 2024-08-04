#!/bin/bash -e

cd ../data

reseek \
  -daliscore_msa 1aab_1j46.afa \
  -input 1aab_1j46.files \
  -output ../output/daliscore_1aab_1j46.tsv \
  -log ../output/daliscore_1aab_1j46.log
