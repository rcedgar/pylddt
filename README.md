# pylddt

Python scripts to calculate DALI Z score and LDDT from a test Multiple Sequence Alignment (MSA) and PDB files with structures. 
Sequences are matched to structures by string matching, the sequence of a PDB file is taken from CA ATOM records.
Structures which do not match the MSA are ignored. Sequences in the MSA which do not match a structure are reported, then ignored.

## Usage 

Scripts are in `src/*.py`. Type `scriptname.py -h` to get usage information.

## To reproduce FoldMason LDDT anomalies

<pre>

cd pylddt/bin
gunzip *.gz

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

<b>Anomalous accounting for gaps with the FoldMason implemention of LDDT</b>.


(a) Three structures used to construct test
cases: (i) scorpion toxin 1sco chain A, (ii) 2iy5 A/a an all-alpha domain at the start of 2iy5 chain A,
and (iii) 4dl1 A/b, an all-beta domain at the start of 4ld1 chain A.


(b) Incorrect alignment of four copies
of 1sco A where the first amino acid is changed to ensure that sequences are distinct, e.g. 1sco/D has
first letter D.

(c) Correct alignment corresponding to (b). 

(d) 1sco A+2iy5 A/a aligned to
1sco A+4ld1 A/b; these structure pairs were manually concatenated to simulate multi-domain proteins
with one common domain. Domains 2iy5 A/a and 4ld1 A/b are not homologous and have no sequence
or structural similarity, therefore (d) is incorrect; these domains should be considered independent
inserts per alignment (e). 

Note that LDDT-fm (FoldMason) assigns 1.0, the highest possible score, to both (b) and (c);
also, LDDT-fm assigns a very low score to correct alignment (e) and a higher score to incorrect
alignment (d). LDDT-mu (Muscle) and Z assign higher scores to the correct alignments.
