#!/bin/bash -e

cd ../test_data

reseek \
  -daliscore_msa 1aab_1j46.afa \
  -input 1aab_1j46.files \
  -output ../test_output/daliscore_1aab_1j46.tsv \
  -log ../test_output/daliscore_1aab_1j46.log
