library(stringr)
library(VennDiagram)
library(ggplot2)
library(dplyr)


import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )%s[^,]|^%s' %s", t_name, t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  # Remove first column and return table.
  t[,c(-1)]
}