# Environmental metagenomics

The workflow and software used in the [Environmental Metagenomics](https://www.genomics-aotearoa.org.nz/projects/environmental-metagenomics) project. Documentation herein does not explain how to install or maintain software. If you are working on the *New Zealand eScience Infrastructure* (NeSI) then most of the software will be already available. Feel free to use the custom scripts and tools here. However, know that they are not actively maintained.

# Assumptions

We assume that you are:
- comfortable with Unix command-line for file system navigation, tools, and syntax
- able to read and understand some bash, Python, and R syntax
- comfortable troubleshooting problems
- understand omics terminology and their synonyms

# Sections

## 1. Sequence processing

Quality check (QC) and preparation of raw reads (DNA and RNA) from Illumina and Oxford Nanopore sequencing platforms from environmental samples or isolates. This is followed by assembly of contigs/scaffolds and assembly QC.

## 2. Prokaryote metagenomics

a. Metagenomic analyses of contigs via gene prediction, annotation, and SSU rRNA sequence reconstruction.

b. Construction, dereplication, and refinement of metagenome-assembled genomes (MAGs) followed by prediction of genomic features such as taxonomic classification, rRNA, tRNA, and CRISPR.

## 3. Viromics

Identification, taxonomic prediction, and functional analyses of:

a. DNA viruses present in metagenomic datasets.

b. RNA viruses present in metatranscriptomic dataset. (Work in progress)

## 4. Transcriptomics

Analyses of metatranscriptomic reads to obtain expression profiles.

---

# Contributors

**Dr David Waite**<br>
The University of Auckland <br>
Genomics Aotearoa<br>
08/10/2019<br>
*Original author of this repository*

**Dr Michael Hoggard**<br>
The University of Auckland<br>
Genomics Aotearoa<br>
06/02/2023<br>
*Viromics and long-read sequencing*

**Jian Sheng Boey**<br>
The University of Auckland<br>
Genomics Aotearoa<br>
*Maintainer of repository*