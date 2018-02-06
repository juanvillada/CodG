################################################################################
################################################################################
###### CodG (C) 2014-2016, Federal University of Viçosa (UFV). All rights reserved.
###### 
######  Z-values is part of CodG
######  
######  Z-values is an R algorithm based on gplots and RColorBrewer packages to 
######  calculate and plot (heatmap) the matrix of squared Z values when determined
######  by other CodG algorithms.
###### 
######## Created by Juan Camilo Villada Arteaga 
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



##################################################################################
################################################################################## 
############################                ###################################### 
############################   Z-values     ######################################
############################                ######################################
################################################################################## 
################################################################################## 

library(gplots)
library(RColorBrewer)

##################################################################################
############################## INPUTS ############################################
##################################################################################
Observed_matrix <- read.csv(file = "Observed_pdCUB.csv")[-1]
mean_matrix <- read.csv(file = "Expected.csv")[-1]
var_matrix <- read.csv(file = "Var.csv")[-1]
codon_rownames <- read.csv("Codon_list_alphabetic.csv", header = F)[[1]]
codon_list <- as.matrix(read.csv("codon_list_Aminoacids_ordered.csv", header = TRUE))
#################################################################################
#################################################################################



################################################################################
############################## FUNCTIONS #######################################
################################################################################
Z_squared <- ((Observed_matrix - mean_matrix)^2) / var_matrix
rownames(Z_squared) <- codon_rownames
j_matrix_squared <- Z_squared
j_matrix_squared <- cbind(Codon = rownames(j_matrix_squared), j_matrix_squared)
j_matrix_squared <- merge(codon_list, j_matrix_squared, by="Codon", all = T)
j_matrix_squared <- j_matrix_squared[order(j_matrix_squared[,2]),]
j_matrix_squared$x <- apply( j_matrix_squared[ , c("Aminoacid", "Codon") ] , 1 , paste , collapse = "-" )
rownames(j_matrix_squared) <- j_matrix_squared[,13]
j_matrix_squared[,1] <- NULL
j_matrix_squared[,1] <- NULL
j_matrix_squared[,11] <- NULL
j_matrix_squared <- as.matrix(j_matrix_squared)
################################################################################
################################################################################


################################################################################
############################### OUTPUTS ########################################
################################################################################
write.csv(j_matrix_squared, "Z-squared_matrix.csv")
#################################################################################
#################################################################################



#################################################################################
#################################### HEATMAP ####################################
#################################################################################

# Palette of orages, used for E. coli heatmap:
my_palette <- colorRampPalette(c("#feedde", "#fdd0a2", "#fdae6b" ,"#fd8d3c", "#f16913", "#d94801", "#8c2d04"))(n = 69)
# Palette of greens, used for yeast plots:
# my_palette <- colorRampPalette(c("#edf8fb", "#ccece6", "#99d8c9" ,"#66c2a4", "#41ae76", "#238b45", "#005824"))(n = 69)

# The subsequent lines of code configure the quadratic plot of the heatmap legend ##########
x <- c(1:7)
y <- (max(j_matrix_squared))/49 * (x^2)
col_breaks <- c(seq(0,y[1],length=10),
                seq((y[1]+0.01),y[2],length=10),
                seq((y[2]+0.01),y[3],length=10),
                seq((y[3]+0.01),y[4],length=10),
                seq((y[4]+0.01),y[5],length=10),
                seq((y[5]+0.01),y[6],length=10),
                seq((y[6]+0.01),y[7],length=10))

# FINAL HEATMAP
heatmap.2(j_matrix_squared,
          density.info="none",
          trace="none",   
          dendrogram="none",
          Rowv = FALSE,
          colsep = c(1:10),
          rowsep = c(1:59),
          Colv=FALSE,
          col = my_palette,
          breaks = col_breaks,
          key = T,
          keysize = 1)
#################################################################################
#################################################################################

