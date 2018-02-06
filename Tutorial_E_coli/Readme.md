# TUTORIAL
## _Escherichia coli_ codon positional dependency

First, you need to configure your working directory in R by the command:
```{r}
setwd("PATH") 
          # Where PATH is the folder containing the genome in Fasta format and all the files of CodG.
``` 
2. Following the schema shown in the Home Page, you just need the coding sequences of the _Escherichia coli_ genome in Fasta format as input. In this case **E_coli_Genome.fasta**. Remember that, because we don't want the high bias in 5'-end regions, and as described in the complete publication, this _E. coli_ genome hasn't the genes coding for proteins with signal peptides.

3. You can start the genomes simulation process by using **SyMuGS** code. The input is the same **E_coli_Genome.fasta**. This step will take a few minutes or hours depending on your system configuration. **SyMuGS** creates the _S_ quantity of simulated genomes, named in this case: _EC_Sim_1.fasta_ , _EC_Sim_2.fasta_ , etc.
You can see an example of this output in the folder: _1.Ecoli_CodG_SyMuGS_

4. Now, by using **QuantiCUB**  you can analyze the original genome (**E_coli_Genome.fasta**) to quantify the codon positional dependency by dividing all the CDSs into ten parts and then counting its codon composition. This algorithm will generate the file _Observed_pdCUB.csv_ which is the observed matrix of codon quantification.
You can see an example of this output in the folder: _2.Ecoli_CodG_QuantiCUB_

5. You can use **ExVar3D** algorithm now to quantify and determine the positional dependency of codons in all the simulated genomes, generating a 3D matrix with the codon quantification of each one of the simulated genomes. Finally, after a couple of munites (or hours), you will have as outputs the files: _Var.csv_ and _Expected.csv_. 
You can see an example of this output in the folder: _3.Ecoli_CodG_ExVar3D_

6. **Z-values** will calculate and save the _Z-squared_matrix.csv_ file and also create the heatmap of the Z-values as follows:
![Zvalues](/Tutorial_E_coli/4.Ecoli_CodG_Z-values/Z-values.png)
You can see an example of this output in the folder: _4.Ecoli_CodG_Z-values_

7.  **Y-values** will calculate and save the _Y-values_matrix.csv_ file and also create the heatmap of the Y-values as follows:
![Yvalues](/Tutorial_E_coli/5.Ecoli_CodG_Y-values/Y-values.png)
You can see an example of this output in the folder: _5.Ecoli_CodG_Y-values_

8. Finally, you can plot the individual codon distribution as a function of the position by using the **IndY-values** algorithm. 
![IndY-values](/Tutorial_E_coli/6.Ecoli_CodG_IndY-values/Codons_Y-value7.png)
You can see an example of this output in the folder: _6.Ecoli_CodG_IndY-values_


9. END! :D
