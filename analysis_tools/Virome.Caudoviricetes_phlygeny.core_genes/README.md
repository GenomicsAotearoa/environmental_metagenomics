# Workflow: Caudoviricetes phylogeny 

Workflow for inference of Caudoviricetes phylogeny via concatenated protein alignments of putative single copy core genes, as used in the study *'Genomic markers show the proliferation of DNA viruses is constrained by adaptation to environmental niche'*.

Method: 

- Based on the method of Low *et al.* (2019) [here](https://doi.org/10.1038/s41564-019-0448-z)
- Previous analyses (Low *et al.* 2019) used 2017 version of VOG and the identified IDs are no longer compatible with the latest version
- in brief, this workflow: 
  - re-identifies putative single copy core genes in Caudoviricetes viruses broadly following the method of Low *et al.* (2019), using latest viralRefSeq references and the latest VOG database
  - extracts the identified putative core genes from viralRefSeq references and Waiwera vOTUs
  - generates filtered concatenated core gene protein alignments for all Caudoviricetes viruses
  - build trees to infer phlyogeny and for visualisation of in iTol

# Software dependencies

- This workflow was developed and tested in Python v3.11.3, and requires the following python libraries: pandas, numpy, re, os, glob, Bio
- Additional software requirements (and tested versions) are outlined in the workflow
- The script compile_dram_annotations.py is available [here](https://github.com/GenomicsAotearoa/environmental_metagenomics/tree/master/3.viromics/a.whole_genome_shotgun/scripts)
