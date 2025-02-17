# Dereplication of putative viral contigs identified by multiple tools

The full workflow for identification of viral contigs from metagenomics assemblies is available [here](https://github.com/GenomicsAotearoa/environmental_metagenomics/blob/master/3.viromics/a.whole_genome_shotgun/1.identify_viral_sequences.md)

Custom scripts used to dereplicate viral contigs identified by multiple tools (VirSorter2; VIBRANT; deepVirFinder), including *summarise_viral_contigs.py* and *virome_per_sample_derep.py*, are available [here](https://github.com/GenomicsAotearoa/environmental_metagenomics/tree/master/3.viromics/a.whole_genome_shotgun/scripts).

# Software dependencies

The custom scripts noted above have been tested with Python v3.8.2, and require the following python libraries: argparse, os, pandas, numpy, re, Bio.SeqIO.FastaIO
