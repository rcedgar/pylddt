for aln in bad_aln correct_aln
do
	../bin/reseek -msta_score $aln.fa -input pdb.files -log reseek.$aln
done

echo
echo
grep LDDT reseek.*
