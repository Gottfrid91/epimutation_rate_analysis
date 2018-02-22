library(data.table)
library(tidyverse)
library(GenomicRanges)
library(LMest)
library(rtracklayer)
library(splitstackshape)
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

#--class distribution--##

#class with lowest fequency
as.data.table(data_dt[context == "CHH", table(L2)])[, min(N)]

#proportions of full data set
full_prop <- as.data.table(data_dt[context == "CHH", table(L2)/dim(data_dt[context=="CHH"])[1]])

#implement the stratified sampling using splitstackshape package
out <- stratified(data_dt[context == "CHH"], c("L2"), 0.001)

#derive proportions
strat_prop <- as.data.table(out[,table(L2)]/dim(out)[1])

#combine the two and see that class proportiant are the same
strat_prop[, full_prop := full_prop[,2]]

#convert into long format
strat_prop_long <- gather(strat_prop, key = "L2")
ggplot(data = strat_prop_long, aes(x=L2, y= value)) + geom_col()
