################################################################################
################################################################################
###### CodG (C) 2014-2016, Federal University of Viçosa (UFV). All rights reserved.
###### 
###### SyMuGS: Synonymous Mutated Genome Simulator
###### 
###### SyMuGS is part of CodG
###### 
###### SyMuGS is an R algorithm based on seqinR library
###### to simulate whole genomes with the insertion of silent mutations. 
###### The algorithm uses a FASTA file of Coding Sequences (CDSs) 
###### as an INPUT, generating the number
###### of genomes indicated by the user. Each simulated genome 
###### will contain random synonymous mutations in all their genes. 
###### It creates a FASTA file of silent mutated CDSs as OUTPUT. 
###### Each simulated genome will have different synonymous mutations, 
###### being different from each other.
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

library(seqinr)

################################################################################
############################## INPUTS ##########################################
################################################################################
# You need to indicate your "Genome.fasta", 
# the file containing the original genome CDSs
real_genome <- read.fasta(file = "Genome.fasta")

# You need to indicate the number of simulated genomes that you want 
# to create, modify the default number (200) for the desired number of genomes.
S <- 200
################################################################################



################################################################################
############################# FUNCTION #########################################
################################################################################

# Function to Simulate and save the new genomes
genomeSimulation <- function(CDSfile, S) {
  for (j in 1:S){
    cat(paste("Wait... Simulating and Saving the Genome # ", j,
              " of ", S, "\n", sep="",collapse=""))
    simulated_genome <- list()
    for (i in 1:length(CDSfile)){
      simulated_genome[i][[1]] <- synsequence(CDSfile[i][[1]])
    }
    write.fasta(simulated_genome, 
                names=getName.list(real_genome), 
                file.out=paste( "SimGenome_", j, ".fasta", 
                                sep = ""))# If you want to specify a directory 
    # to save your simulated genomes
    # change the "SimGenome_" for 
    # "NEW_PATH/SimGenome_" 
  }
  cat(paste("***** FINISHED!", "*****\n ||| ", S,
            " ||| Genomes Simulated and Saved", sep="",collapse=""))
}
################################################################################
################################################################################


################################################################################
############################### OUTPUTS ########################################
################################################################################
########## Execute the following line to start the genome simulations ##########
genomeSimulation(real_genome, S)
################################################################################
################################################################################


