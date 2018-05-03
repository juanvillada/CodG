################################################################################
###### CodG (C) 2014-2016, Federal University of Viçosa (UFV). All rights reserved.
###### 
######  IndY-values is part of CodG
######  
######  IndY-values is an R algorithm to 
######  calculate and plot the individual codon distribution of Y-values 
######  among different positions of the genes divided in Bins
######  when determinedby other CodG algorithms.
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


# IndY-values 
# INDIVIDUAL CODON DISTRIBUTION OF Y-VALUES THROUGH BINS

# INPUT
j_matrix <- read.csv("Y-values_matrix.csv", row.names = 1)


# FUNCTIONS & PDF OUTPUTS
b <- 1
c <- 4
for (n in 1:15){
  pdf(file = paste("Codons_Y-value",n,".pdf"), useDingbats=FALSE)
  par(mfrow=c(2,2))
  for (u in b:c){
    x <- c(1:10)
    y <- j_matrix[u,]
    plot(c(1:10), j_matrix[u,], ylim = c((min(j_matrix)), (max(j_matrix))), pch=16, xlab = "Bin", ylab = "Y-value",
         main = rownames(j_matrix)[u])
    lines(loess.smooth(x, y, col="#bdbdbd", lty=1, lwd=4, family = "gaussian"))
    #lines(x, y, col="#bdbdbd", lty=1, lwd=4)
    abline(h=0, lwd=1, lty=2, col="black")
  }
  dev.off()
  b <- b+4
  if (c == 56) {
    c <- c+3
  } else {
    c <- c+4
  }
}
