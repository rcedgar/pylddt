foldmason=/z/sw/foldmason/foldmason/bin/foldmason

rm -rf tmp/

cd ../foldmason

$foldmason easy-msa ../pdbdir_pp/ pp tmp/
$foldmason easy-msa ../pdbdir_fl/ fl tmp/

rm -rf tmp/

/bin/cp -v pp_aa.fa ../msa_pp/foldmason
/bin/cp -v fl_aa.fa ../msa_fl/foldmason
