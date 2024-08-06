for pp in pp fl
do
	../../src/lddt_foldmason.py --msa ../msa_$pp/foldmason --pdbfiles ../files/$pp.files 
	../../src/lddt_foldmason.py --msa ../msa_$pp/muscle_mega --pdbfiles ../files/$pp.files 
	../../src/lddt_foldmason.py --msa ../msa_$pp/muscle_aa --pdbfiles ../files/$pp.files 

	../../src/daliz.py --msa ../msa_$pp/foldmason --pdbfiles ../files/$pp.files 
	../../src/daliz.py --msa ../msa_$pp/muscle_mega --pdbfiles ../files/$pp.files 
	../../src/daliz.py --msa ../msa_$pp/muscle_aa  --pdbfiles ../files/$pp.files 
done | tee ../scores/scores.txt
