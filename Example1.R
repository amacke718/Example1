# This is Michal's Example
# Other Notes:
# "c" means "concatenate", which commonly means to link/chain objects together, which is what you do when you use c()

file <- "TCGA_breast_cancer_LumA_vs_Basal_PAM50.tsv"
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold <- 3

header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

# dim(data) shows the dimension of the data 49 Columns & 329 Rows, but really thats including the headers, so it is 328 Rows
# data$sample_id=='name of observation' will allow you to filter your data frame - show a particular header/value 
# data[data$sample_id=='BCL2',] will pull out every data point in accordance with your filter
# The use of the comma and space here, indicates "ROW FILTERING"
# View(data) {The capital V here matters} - opens the data in a window
# data[,c(6)]
LumA <- data[data$sample_id %in% first10,header2=='Luminal A']
Basal <- data[data$sample_id %in% first10,header2=='Basal-like']

LumA_groups <- split(colnames(LumA), sample(1:nfold, ncol(LumA), replace=T))
Basal_groups <- split(colnames(Basal), sample(1:nfold, ncol(Basal), replace=T))
  
result <- array()
  
for (test_group in 1:nfold) {
    
    testA <- LumA[,colnames(LumA) %in% unlist(LumA_groups[test_group])]
    testB <- Basal[,colnames(Basal) %in% unlist(Basal_groups[test_group])]
    
    trainingA <- LumA[,!(colnames(LumA) %in% unlist(LumA_groups[test_group]))]
    trainingB <- Basal[,!(colnames(Basal) %in% unlist(Basal_groups[test_group]))]
    
    centroidA <- rowMeans(trainingA)
    centroidB <- rowMeans(trainingB)
    
    misclassifiedA <- sum(sapply(testA, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))>0 }))
    misclassifiedB <- sum(sapply(testB, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))<0 }))
    
    result[test_group] <- (misclassifiedA+misclassifiedB)/(ncol(testA)+ncol(testB))
}
  
mean(result)
sd(result)

