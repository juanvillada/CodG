################################################################################
###### CodG (C) 2014-2016, Federal University of Viçosa (UFV). All rights reserved.
###### 
######  Y-values is part of CodG
######  
######  Y-values is an R algorithm based on gplots and RColorBrewer packages to 
######  calculate and plot (heatmap) the matrix of Y-values when determined
######  by other CodG algorithms.
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

# Y-values

library(gplots)
library(RColorBrewer)

# INPUTS
Observed_matrix <- read.csv(file = "Observed_pdCUB.csv")[-1]
mean_matrix <- read.csv(file = "Expected.csv")[-1]
codon_rownames <- read.csv("Codon_list_alphabetic.csv", header = F)[[1]]
codon_list <- as.matrix(read.csv("codon_list_Aminoacids_ordered.csv", header = TRUE))

# FUNCTIONS
Y <- log2(Observed_matrix/mean_matrix)
rownames(Y) <- codon_rownames
j_matrix <- Y
j_matrix <- cbind(Codon = rownames(j_matrix), j_matrix)
j_matrix <- merge(codon_list, j_matrix, by="Codon", all = T)
j_matrix <- j_matrix[order(j_matrix[,2]),]
j_matrix$x <- apply( j_matrix[ , c("Aminoacid", "Codon") ] , 1 , paste , collapse = "-" )
rownames(j_matrix) <- j_matrix[,13]
j_matrix[,1] <- NULL
j_matrix[,1] <- NULL
j_matrix[,11] <- NULL
j_matrix <- as.matrix(j_matrix)


# OUTPUTS
write.csv(j_matrix, "Y-values_matrix.csv")

# HEATMAP
my_palette2 <- colorRampPalette(c("#404040", "#bababa", "#ffffff", "#f4a582", "#ca0020"))(n = 99)
heatmap.2(j_matrix,
          density.info="none",  # turns off density plot inside color legend
          trace="none",   # enable color transition at specified limits
          dendrogram="none",
          Rowv = FALSE,
          sepcolor = "white",
          colsep = c(1:10),
          rowsep = c(1:59),
          col = my_palette2,
          Colv=FALSE)
