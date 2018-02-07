# CodG: Analyzing codon positional dependency from genome-scale data

https://doi.org/10.1093/dnares/dsx014

![CodG Schema](/Images/CodG.png)

____

# Click [here](https://github.com/juanvillada/CodG/tree/master/Tutorial_E_coli) for the tutorial.

____

# CodG

## Description
`CodG` is a set of scripts developed in R language that automate the analysis of preferred codon usage as a function of the position within genes using genomic data.

This program is based on the `SeqinR` and `BioStrings` packages. The `CodG` has three main scripts: `SyMuGS` (Synonymous Mutated Genome Simulator), `QuantiCUB` (Position-dependent Codon Usage Bias Quantifier), `ExVar3D` (Expected Value and Variance calculator); and other scripts for support visualization and other analyses: `Z-values`, `Y-values` and `IndY-values`.

### SyMuGS: Synonymous Mutated Genome Simulator
The algorithm, `SyMuGS`, simulates complete genomes introducing random synonymous mutations in each of the coding sequences. The preferred codon usage bias of each of the sequences is conserved. This algorithm uses the file containing coding sequences in FASTA format as input and generates as output a different file in FASTA format for each simulated genome.

### QuantiCUB: Position-dependent Codon Usage Bias Quantifier
`QuantiCUB` quantifies codon usage bias as a function of the position according to the length of each coding sequence in the genome, using as input the FASTA file format containing the coding sequences of the genome to be analyzed. It generates the _Observed_pdCUB_ array of 10 positions for 59 codons in a comma separated values ​​(CSV) format file.

### ExVar3D: Expected Value and Variance calculator
The `ExVar3D` algorithm calculates the position-dependent codon usage bias for each of the multiple simulated genomes, having as input the sequences of the simulated genomes in FASTA format which are obtained with the `SyMuGS` script. The program is based on the `QuantiCUB` script, but other than this, the `ExVar3D` generates a virtual tridimensional matrix with the quantification of codons of each genome, and from this generates another two arrays: one is an expected value of the matrix of each codon by position, and an array of variance codons for each position. The output format of these two matrices is CSV.

### Z-values, Y-values and IndY-values
`Z-values`, `Y-values` and `IndY-values` are algorithms implemented to create comprehensive figures in order to illustrate the results obtained from the previous analyses.

`CodG` compares the codon distribution differences of the original genome with a large number of simulated genomes, allowing for the analysis of the evolutionary selective forces acting on individual codons. This set of scripts can be useful for the analysis of microbial genes and genomes using different approaches such as in evolutionary research, expression of recombinant proteins and synthetic biology.

____

_Please contact [Juan C. Villada](mailto:juan.arteaga@ufv.br) if you have any questions_
____