# pylddt

Python scripts to calculate DALI Z score and LDDT from a test Multiple Sequence Alignment (MSA) and PDB files with structures. 
Sequences are matched to structures by string matching, the sequence of a PDB file is taken from CA ATOM records.
Structures which do not match the MSA are ignored. Sequences in the MSA which do not match a structure are reported, then ignored.

## Usage 

Scripts are in `src/*.py`. Type `scriptname.py -h` to get usage information.

## To reproduce FoldMason LDDT anomalies

<pre>
cd pylddt/scorpion_toxin
./run_msa2lddt.bash
./run_reseek.bash

cd pylddt/scorpion_toxin3
./run_msa2lddt.bash
./run_reseek.bash
</pre>

## References

Mariani, Valerio, et al. "LDDT: a local superposition-free score for comparing protein structures and models using distance difference tests." <i>Bioinformatics</i> 29.21 (2013): 2722-2728.

Holm, Liisa, and Chris Sander. "Protein structure comparison by alignment of distance matrices." <i>Journal of molecular biology</i> 233.1 (1993): 123-138.

Gilchrist CL, Mirdita M, Steinegger M. Multiple Protein Structure Alignment at Scale with FoldMason. bioRxiv. 2024 Aug 1:2024-08.


![FoldMason LDDT anomalies](https://github.com/rcedgar/pylddt/raw/main/results/FoldMason_LDDT_anomaly_figure.png)
