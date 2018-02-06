################################################################################
################################################################################
###### CodG (C) 2014-2016, Federal University of Viçosa (UFV). All rights reserved.
###### 
######  QuantiCUB: Position-dependent Codon Usage Bias Quantifier
######  
######  QuantiCUB is part of CodG
######  
######  QuantiCUB is an R algorithm based on BioStrings package to quantify 
######  the position-dependant Codon Usabe Bias relative to each 
######  Coding Sequence (CDS) length. The algorithm uses a FASTA file 
######  of CDSs as an INPUT. CDSs with length lower than 120 nucleotides 
######  are excluded. Each CDS is partitioned in 10 parts (bins), 
######  the start and stop codons are excluded and each codon is quantified
######  by gene, then it is saved the total quantification of all the gene 
######  codons being part of the same bin. Thus, it is generated a matrix 
######  (59x10) of 59 codons (the non-redundant codons are excluded) 
######  and 10 bins. 
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
################################################################################

library(Biostrings)

################################################################################
############################## FUNCTIONS #######################################
################################################################################

relative_pD <- function(CDS_file, Q){
  CDS_file <- subseq(CDS_file, start=4, end = -4) 	# Start and Stop 
  													# codon exclusion
  													
  CDS_file <- CDS_file[width(CDS_file)>=120] 	# CDSs with length lower 
  												# than 120 are excluded
  												
  W <- round(((width(CDS_file))/3)/Q)*3 # Vector of each CDS bin size
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
  P <- P[-c(15,49,51,57,59),]   # Removing non redundant 
  								# codons (ATG,TAA,TAG,TGA,TGG)
  return(P)
  }
  

################################################################################
############################## INPUTS ##########################################
################################################################################

CDS_file <- readDNAStringSet(filepath = "Genome.fasta") 	# Indicate the 
# PATH of the 
# CDS FASTA file

Q <- 10 #Number of relative bins
################################################################################



################################################################################
########################## START THE PROCESS ###################################
################################################################################
Observed_matrix <- relative_pD(CDS_file = CDS_file, Q)
################################################################################



################################################################################
############################### OUTPUT #########################################
################################################################################
write.csv(Observed_matrix, file = "Observed_pdCUB.csv") # Indicate the 
														# PATH to save 
														# the matrix of 
														# Codon 
														# Quantification	
################################################################################
################################################################################
