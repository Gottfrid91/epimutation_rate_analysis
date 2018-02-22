library(data.table)
library(tidyverse)
library(GenomicRanges)
library(splitstackshape)

handle_data <- function(data_dt, cont){
  '
  Input: 
  data_dt: Full data frame containing all data
  context: context of cytocins 
  Output: 
  S: Matrix containing data ready for EM algorithm and Markov chain fitting   
  '
  
  #draw a sub-sample of observations
  cg_sub_dt <- stratified(data_dt[context == cont], c("L2"), 0.001)
  n <- nrow(cg_sub_dt)
  
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
  
  return(S)
}

#number to binary function 
number2binary = function(number) {
  binary_vector = rev(as.numeric(intToBits(number)))
  paste0(binary_vector[-(1:(length(binary_vector) - 6))], collapse = '')
}

# Function to calculate first-order Markov transition matrix.
# Each *row* corresponds to a single run of the Markov chain
trans.matrix <- function(X, prob=T)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  if(prob) tt <- tt / rowSums(tt)
  tt
}

