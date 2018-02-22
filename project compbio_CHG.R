library(data.table)
library(tidyverse)
library(GenomicRanges)
library(splitstackshape)
library(Publish) 
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
cont = "CHG"

#get data into preprocessed and matrix form
S <- handle_data(data_dt, cont)

##-- Start deriving arguments for EM algortithm --##

#one cycle transition counts
N_1 <- trans.matrix(S, prob=F)
#two cycle transition matrix
S_2 <- S[,c("V2","V4")]
N_2 <- trans.matrix(S_2, prob=F)
#three cycle transition matrix
S_3 <- S[,c("V5","V8")]
N_3 <- trans.matrix(S_3, prob=F)

#three different initial transition matrix 
M_1 <- trans.matrix(S, prob=T)
M_2 <- trans.matrix(S_2, prob=T)**(0.5)
M_3 <- trans.matrix(S_3, prob=T)**(0.33)

#threshold for EM
eps <- 0.01
T_old <- matrix(c(0,0,0,0), ncol = 2, nrow = 2)

em_list_CHG <- EMalgorithm(M_1, N_1, N_2, N_3, eps)
T_2_CHG <- EMalgorithm(M_2, N_1, N_2, N_3, eps)
T_3_CHG <- EMalgorithm(M_3, N_1, N_2, N_3, eps)

#call bootstrap
result_CHG <- bootstrap(em_list_CHG, N_1, N_2, N_3, 1000)
#save working data
saveRDS(result_CHG, file = 'result_CHG.RData')
