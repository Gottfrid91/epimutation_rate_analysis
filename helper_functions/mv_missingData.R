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

data_2 <- load(data_file_context_2)

data_3 <- readRDS(data_file_context_3)
data_3

import.gff3()

number2binary = function(number) {
  binary_vector = rev(as.numeric(intToBits(number)))
  paste0(binary_vector[-(1:(length(binary_vector) - 6))], collapse = '')
}

data_dt <- data.table(context = data_3@elementMetadata[,1], 
                      L2 = data_3@elementMetadata[,2], 
                      L8 = data_3@elementMetadata[,3])

#draw a random sub-sample of 10000 observations
n <- 10000

markovmodel <- function(s, n){
  #get random sub_sample
  cg_sub_dt <- sample_n(data_dt[context == context], n)
  
  #
  string <- list()
  k <- 1
  for(i in 1:n){
    string_i <- number2binary(cg_sub_dt[i,2])
    string[[k]] <- string_i
    k <- k+1
  }
  
  
  cg_sub_dt[,methylation :=unlist(string)]
  
  cg_sub_dt[  by = unique(methylation)]
  #loop to split string of generational methylations into own columns
  length(string)
  dt <- matrix(0,n,6)
  
  for(i in 1:n){
    dt[i,] <- rbind(strsplit(cg_sub_dt$methylation,"")[[i]])
  }
  
  #make entries numeric
  dt <- matrix(sapply(dt, as.numeric), n, 6)
  
  #creating markov model data for a context CG
  mod <- 1
  k <- 2
  start <- 1
  
  out <- aggr_data(dt)
  yv <- out$freq
  S <- out$data_dis
  
  S <- as.data.table(S)
  V_i = rep(NA, dim(S)[1])
  S = as.matrix(data.table(V0 = S[,V1], V1 = V_i, V2 = S[,V2], V3 = S[,V3], V4 = S[,V4], 
                           V5 = S[,V5], V6 = V_i, V7 = V_i, V8 = S[,V6]))
  
  mc <- est_lm_basic(S = S,yv = yv, k = 2, start = 0, mod = 1,out_se = TRUE)
  
  return(mc)
}
n <- 10000
context <- "CG"
mc_cg <- markovmodel(context, n)
context <- "CHH"
mc_chh <- markovmodel(context, n)
context <- "CHG"
mc_chg <- markovmodel(context, n)

#Transition probabilities:
T_cg <- mc_cg$Pi[,,6]
T_chh <- mc_chh$Pi[,,6]
T_chg <- mc_chg$Pi[,,6]
#T_all <- mc_all$Pi[,,6]

#Standard errors for transition probabilities:
stderr_cg <- mc_cg$sePi[,,6]
stderr_chh <- mc_chh$sePi[,,6]
stderr_chg <- mc_chg$sePi[,,6]
#stderr_all <- mc_all$sePi[,,6]