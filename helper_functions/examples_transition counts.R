library(data.table)
library(tidyverse)
library(GenomicRanges)
library(splitstackshape)
source('/WORKING_DIR_HERE/stat_functions.R')
source('/WORKING_DIR_HERE/data_helper_functions.R')

#get full data
DATA_DIR_3 <- 'WORKING_DIR_HERE/data/'
data_file_context_3 <- file.path(DATA_DIR_3, 'methpatterns_whole_genome.RData')
#read S4 data object
data_3 <- readRDS(data_file_context_3)

#extract wanted column from full data set and convert to data table
data_dt <- data.table(context = data_3@elementMetadata[,1], L2 = data_3@elementMetadata[,2],
                      L8 = data_3@elementMetadata[,3])

#set context and call prerocessing function
cont = "CHH"

#get data into preprocessed and matrix form
S <- handle_data(data_dt, cont)

trans.matrix <- function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}


##--Verify that counting works --##

X_2 <- head(S,2)
print(X_2)
#verify the count
trans.matrix(X_2, prob = F)

#another row added to matrix
X_3 <- head(S,3)
print(X_3)
trans.matrix(X_3, prob = F)

#verify counts when it does not only contain 0's
X_4 <- head(S[S[, "V1"] == 1 & S[, "V0"] == 1,])
print(X_4)
trans.matrix(X_4, prob = F)


