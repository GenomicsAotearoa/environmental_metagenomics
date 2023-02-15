# Metagenomic sequence analyses

## Introduction

After assembling your reads, we can attempt to characterise "what" your samples are capable of doing (functional gene annotation), and "who" is in them (SSU rRNA reconstruction and taxonomic classification). 

To perform the former, we need to
1. predict (full or partial) protein-coding genes 
2. annotate the predictions against large sequence or protein domain databases via sequence alignment or Hidden Markov Models (HMM)

To perform the latter, we need to
1. reconstruct a full-length SSU rRNA sequence from quality trimmed reads
2. classify the reconstructed sequence against the SILVA database

Along with appropriate read mapping, the resulting table from SSU rRNA reconstruction can be interpreted and analysed like an OTU/ASV table.

## Required data

- Quality trimmed reads
- Assembled contigs/scaffolds
- Downloaded databases
