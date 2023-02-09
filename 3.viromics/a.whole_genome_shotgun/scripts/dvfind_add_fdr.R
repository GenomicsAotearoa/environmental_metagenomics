#!/usr/bin/env Rscript

# Load arguments
args = commandArgs(trailingOnly=TRUE)

# Load libraries
library(readr)
library(dplyr)
library(qvalue)

# Load file, calculate qvalues, calculate fdr adjusted p-values, sort by p.adj
## Filter with parameters: score >= 0.9 & pvalue <= 0.05 & p.adj <= 0.1
## export as tsv (overwrite original)
result <- read_tsv(args[1]) %>%
    mutate(., qvalue = tryCatch(qvalue(.$pvalue, pi0.meth="bootstrap")$qvalues, error=function(e) "error")) %>%
    mutate(., p.adj = p.adjust(.$pvalue, method="fdr")) %>%
    arrange(., p.adj) %>%
    filter(., (score >= 0.9 & pvalue <=0.05 & p.adj <= 0.1)) %>%
    write_tsv(., args[1])
