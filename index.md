

![CodG Schema](/Images/CodG.png)

____



### _Please note that the National Institute for Industrial Property ([INPI](http://www.inpi.gov.br/)) of the Brazilian Government has determined the rights of distribution of this software (Process: BR512016000508-4). Please contact [Juan Villada](mailto:juan.arteaga@ufv.br) if you want to use this software for academic research. Please do your requests using your institutional e-mail._

____

# CodG

## Description
CodG is a program developed in R language that automates the testing procedures of the preferred codon usage as a function of the position within genes using genome-scale data. 

This program is based on the SeqinR and BioStrings packages. The program uses six scripts: SyMuGS (Synonymous Mutated Genome Simulator), QuantiCUB (Position-dependent Codon Usage Bias Quantifier), ExVar3D (Expected Value and Variance calculator), Z-values, Y-values and IndY-values.

### SyMuGS: Synonymous Mutated Genomes Simulator
The script, SyMuGS, simulates complete genomes introducing random synonymous mutations in each of the coding sequences. The preferred codon usage bias of each of the sequences is conserved. This script uses the file containing coding sequences in FASTA format as input and generates as output, a different file in FASTA format for each simulated genome.

### QuantiCUB: Position-dependent Codon Usage Bias Quantifier
QuantiCUB quantifies codon usage bias as a function of the position according to the length of each coding sequence in the genome, using as input the FASTA file format containing the coding sequences of the genome to be analyzed. It generates the “Observed_pdCUB” array of 10 positions for 59 codons in a comma separated values (CSV) format.

### ExVar3D: Expected Value and Variance calculator
The ExVar3D script calculates the position-dependent codon usage bias for each of the multiple simulated genomes, having as input the sequences of the simulated genomes in FASTA format which are obtained with the SyMuGS script. The program is based on the QuantiCUB script, but other than this, the ExVar3D generates a virtual tridimensional matrix with the quantification of codons of each genome, and from this generates another two arrays: one is an expected value of the matrix of each codon by position, and an array of variance codons for each position. The output format of these two matrices is CSV.

### Z-values, Y-values and IndY-values
Z-values, Y-values and IndY-values are scripts implemented to create figures to illustrate the results obtained from the  previous analyses.

CodG compares the codon distribution differences of the original genome with a large number of simulated genomes, allowing for the analysis of the evolutionary selective forces acting on individual codons. This suite of scripts is a useful tool for the analysis of genes and genomes from microorganisms using different approaches such as evolutionary research, expression of recombinant proteins and synthetic biology.
