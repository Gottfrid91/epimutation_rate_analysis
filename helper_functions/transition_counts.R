library(data.table)
library(tidyverse)
library(GenomicRanges)
library(LMest)
library(rtracklayer)

source('/WORKING_DIR_HERE/stat_functions.R')
source('/WORKING_DIR_HERE/data_helper_functions.R')

#get full data
DATA_DIR_3 <- 'WORKING_DIR_HERE/data/'
data_file_context_3 <- file.path(DATA_DIR_3, 'methpatterns_whole_genome.RData')

#read in S4 data object
data_3 <- readRDS(data_file_context_3)

#number to bi nary function 
number2binary = function(number) {
  binary_vector = rev(as.numeric(intToBits(number)))
  paste0(binary_vector[-(1:(length(binary_vector) - 6))], collapse = '')
}

data_dt <- data.table(context = data_3@elementMetadata[,1], L2 = data_3@elementMetadata[,2],
                      L8 = data_3@elementMetadata[,3])

#draw a random sub-sample of 10000 observations
n <- 10000
cg_sub_dt <- sample_n(data_dt[context == "CG"], n)

#convert number to binary and make into string
string <- list()
k <- 1
for(i in 1:n){
  string_i <- number2binary(cg_sub_dt[i,2])
  string[[k]] <- string_i
  k <- k+1
}

#unlist string to make into a matrix
cg_sub_dt[,methylation :=unlist(string)]

#loop to split string of generational methylations into own columns
length(string)

#initialize matrix object to fill
dt <- matrix(0,n,6)

#split string and insert in to dt matrix object
for(i in 1:n){
  dt[i,] <- rbind(strsplit(cg_sub_dt$methylation,"")[[i]])
}

#make entries numeric
dt <- matrix(sapply(dt, as.numeric), n, 6)

#insert NA's
S <- as.data.table(dt)
V_i <- rep(NA, dim(S)[1])
S <- as.matrix(data.table(V0 = S[,V1], V1 = S[,V2], V2 = S[,V3], V3 = V_i, V4 = S[,V4], 
                         V5 = S[,V5], V6 = V_i, V7 = V_i, V8 = S[,V6]))

# Function to calculate first-order Markov transition matrix.
# Each *row* corresponds to a single run of the Markov chain
trans.matrix <- function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}

#one cycle transition counts
N_1 <- trans.matrix(S, prob=F)

#two cycle transition matrix
S_2 <- S[,c("V2","V4")]
N_2 <- trans.matrix(S_2, prob=F)

#three cykle transition matrix
S_3 <- S[,c("V5","V8")]
N_3 <- trans.matrix(S_3, prob=F)

