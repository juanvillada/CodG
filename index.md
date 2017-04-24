

## _Please note that the National Institute for Industrial Property ([INPI](http://www.inpi.gov.br/)) of the Brasilian Goverment has determined the rights of distribution of this software. Please contact [Juan Villada](juan.arteaga@ufv.br) if you want to use the first release of this software for academic research. You can request a copy of the scripts using your institutional e-mail._

____

![CodG Schema](/Images/CodG.png)

____

# CodG

## Description
CodG is a program developed in R language that automates the testing procedures of the preferred codon usage as a function of the position within genes using genome-scale data. 

This program is based on the SeqinR and BioStrings packages. The program uses six algorithms: SyMuGS (Synonymous Mutated Genome Simulator), QuantiCUB (Position-dependent Codon Usage Bias Quantifier), ExVar3D (Expected Value and Variance calculator), Z-values, Y-values and IndY-values.

### SyMuGS: Synonymous Mutated Genome Simulator
The algorithm, SyMuGS, simulates complete genomes introducing random synonymous mutations in each of the coding sequences. The preferred codon usage bias of each of the sequences is conserved. This algorithm uses the file containing coding sequences in FASTA format as input and generates as output, a different file in FASTA format for each simulated genome.

### QuantiCUB: Position-dependent Codon Usage Bias Quantifier
QuantiCUB quantifies codon usage bias as a function of the position according to the length of each coding sequence in the genome, using as input the FASTA file format containing the coding sequences of the genome to be analyzed. It generates the “Observed_pdCUB” array of 10 positions for 59 codons in a comma separated values ​​(CSV) format.

### ExVar3D: Expected Value and Variance calculator
The ExVar3D algorithm calculates the position-dependent codon usage bias for each of the multiple simulated genomes, having as input the sequences of the simulated genomes in FASTA format which are obtained with the SyMuGS algorithm. The program is based on the QuantiCUB algorithm, but other than this, the ExVar3D generates a virtual tridimensional matrix with the quantification of codons of each genome, and from this generates another two arrays: one is an expected value of the matrix of each codon by position, and an array of variance codons for each position. The output format of these two matrices is CSV.

### Z-values, Y-values and IndY-values
Z-values, Y-values and IndY-values are algorithms implemented to create comprehensive figures in order to illustrate the complex results obtained of the entire previous experimentation.

Thus, CodG compares the codon distribution differences of the original genome with a large number of simulated genomes, allowing for the analysis of the evolutionary selective forces acting on individual codons. This program is a useful tool for the analysis of genes and genomes from microorganisms using different approaches such as evolutionary research, expression of recombinant proteins and synthetic biology.
