################################################################################
###### CodG (C) 2014-2016, Federal University of Viçosa (UFV). All rights reserved.
###### 
######  ExVar3D: Expected Value and Variance calculator
######  
######  ExVar3D is part of CodG
######  
######  ExVar3D is an R algorithm based on SeqinR and BioStrings packages to 
######  quantify position-dependent codon usage bias relative to each 
######  CDS length as described for QuantiCUB, but implemented for a 
######  high number of simulated genomes. This algorithm allows the user
######  to contruct a 3D matrix of dimensions 59 (Codons), 10 (Bins), 
######  G (number of simulated Genomes). Then, the algorithm calculates 
######  the Expected value and Variance for each codons by bin (position). 
######  The user needs to indicate the PATH containing the Simulated 
######  Genomes, the quantity of Simulated Genomes as an INPUT and then 
######  it is created two matrixes as OUTPUT with the Expected value and 
######  Variance of the codon quantification related to the genomes with 
######  silent mutaions in their CDSs.
###### 
####### Created by Juan Camilo Villada Arteaga 
####### supervised by Wendel Batista da Silveira 
####### Department of Microbiology 
####### Federal University of Viçosa (UFV) - Brazil. 2016.
####### 
####### CodG has been deposited under the process number BR 51 2016 000508 4
####### into The National Institute of Industrial Property, 
####### more commonly abbreviated to INPI, is the intellectual property office 
####### responsible for industrial property in Brazil. 
####### 
####### The UFV hereby grants to you the non-exclusive right to use the CodG 
####### software only, solely for non-commercial, research-only purposes.
################################################################################

library("seqinr")
library("Biostrings")


# FUNCTIONS
relative_pD <- function(CDS_file, Q){
  CDS_file <- subseq(CDS_file, start=4, end = -4)
  CDS_file <- CDS_file[width(CDS_file)>=120]
  W <- round(((width(CDS_file))/3)/Q)*3
  Y <- rep(0,64)
  X <- vector()
  P <- matrix(nrow = 64, ncol = Q)
  j <-1
  for (k in 0:(Q-2)){
    for (i in 1:(length((CDS_file)))){
      X <- trinucleotideFrequency(CDS_file[[i]][((k*W[i])+1):
      								((k+1)*W[i])],step = 3)
      Y <- X+Y
    }
    P[,j] <- Y
    Y <- rep(0,64)
    j <- j+1
  }
  for (s in 1:(length(CDS_file))) {
    X <- trinucleotideFrequency(CDS_file[[s]][(((Q-1)*W[s])+1):
    								(width(CDS_file[s]))],step = 3)
    Y <- X+Y
  }
  P[,Q] <- Y
  P <- P[-c(15,49,51,57,59),]   #Removing non redundant codons (ATG,TAA,TAG,TGA,TGG)
  return(P)
}

# INPUTS
G <- 200 # Indicate Quantity of Simulated Genomes to analyse
Q <- 10 # Bins

# CREATIING THE 3D MATRIX 
m3D <- array(NA, c(59,Q,G))
for (k in 1:G){

  cat(paste("Wait... analysing the Genome # ", k, " of ", 
  			G, "\n", sep="",collapse=""))
  			
  CDS_file <- readDNAStringSet(filepath = paste("EC_Sim_", k, ".fasta", sep = "")) #Indicate in "SimGenome_" the PATH (PATH/"SimGenome_") of the simulated genomes created by SyMuGS
  								
  M <- relative_pD(CDS_file = CDS_file, Q)
  
  m3D[,,k] <- M # Contruction of the 3D matrix of codon quantification as used in QuantiCUB 
}

# CREATING THE EXPECTED(MEAN) AND VARIANCE MATRIXES
# Construction of the expected matrix of the total of quantified genomes by bin and codon
mean_matrix <- matrix(nrow = 59, ncol = Q)

# Construction of the variance value matrix of the total of quantified genomes by bin and codon
var_matrix <- matrix(nrow = 59, ncol = Q)

										
gg <- vector()
for (j in 1:Q){
  for (i in 1:59){
    for (k in 1:G){
      gg[k] <- c(m3D[i,j,k])
      mean_matrix[i,j] <- mean(gg)
      var_matrix[i,j] <- var(gg)
    }
  }
}


# OUTPUTS
# Indicate the PATH to save the Expected Values Matrix
write.csv(mean_matrix, file = "Expected.csv")


# Indicate the PATH to save the Variance Values Matrix
write.csv(var_matrix, file = "Var.csv")
