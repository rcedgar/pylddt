# pylddt

Python script to calculate DALI Z score and LDDT from a test Multiple Sequence Alignment (MSA) and PDB files with structures.
Reports Z and LDDT for each pair of sequences, and the mean Z and LDDT over all pairs. Sequences are matched to structures by
string matching, the sequence of a PDB file is taken from CA ATOM records.Structures which do not match the MSA are
ignored. Sequences in the MSA which do not match a structure are reported, then ignored.

## Usage 
<pre>
pylddt.py [-h] --msa MSA --pdbfiles PDBFILES [--radius RADIUS] [--dists DISTS] [--horizon HORIZON]
                 [--diagwt DIAGWT] [--output OUTPUT]

  -h, --help           show this help message and exit
  --msa MSA            Multiple Sequence Alignment (FASTA format)
  --pdbfiles PDBFILES  Text file with one PDB pathname per line
  --radius RADIUS      LDDT inclusion radius (default 15)
  --dists DISTS        LDDT distance thresholds, comma-separated (default 0.5,1,2,4)
  --horizon HORIZON    DALI horizon (default 20)
  --diagwt DIAGWT      DALI diagonal weight (default 0.2)
</pre>

## References

Mariani, Valerio, et al. "lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests." <i>Bioinformatics</i> 29.21 (2013): 2722-2728.

Holm, Liisa, and Chris Sander. "Protein structure comparison by alignment of distance matrices." <i>Journal of molecular biology</i> 233.1 (1993): 123-138.
