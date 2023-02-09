#!/usr/bin/env python

'''virome_per_sample_derep.py

Dereplicate putative viral contigs identified by multiple tools and generate per-sample dereplicated fasta file 

- Requires assembly fasta file and summary_table_VIRUSES.txt as inputs 
- If prophage identified (and excised) via VIBRANT or VirSorter2, the output fasta files from these tools are required to extract the prophage sequences for the dereplicated output
- Script also filters out any included full contigs that were also identified as (excised) prophage via VIBRANT or VirSorter2 (i.e. only retaining trimmed/excised versions)


Required parameters:

--assembly_fasta (-a) : Assembly fasta file used as input for VIBRANT, VirSorter2, or DeepVirFinder
--summary_table (-t) : summary_table_VIRUSES.txt (generated via wgs_summarise_viral_contigs.py)

Optional parameters:

--vibrant (-v) : Fasta file of putative phage output from vibrant: SampleX.phages_combined.fna
--virsorter2 (-s) : Fasta file of putative phage output from VirSorter2: SampleX-final-viral-combined.fa
--output (-o) : Dereplicated fasta filename. Default = 'per_sample_derep.fasta'.
'''

from argparse import ArgumentParser
import os
import pandas as pd
import numpy as np
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument('-a', '--assembly_fasta', dest="assembly_fasta",
                    help="Assembly fasta file used as input for VIBRANT, VirSorter2, or DeepVirFinder", required=True)
parser.add_argument('-t', '--summary_table', dest="summary_table",
                    help="summary_table_VIRUSES.txt (generated via wgs_summarise_viral_contigs.py)", required=True)
parser.add_argument('-v', '--vibrant', dest="vibrant_fasta",
                    help="Fasta file of putative phage output from vibrant: SampleX.phages_combined.fna", default=None)
parser.add_argument('-s', '--virsorter2', dest="vsort2_fasta",
                    help="Fasta file of putative phage output from VirSorter2: SampleX-final-viral-combined.fa", default=None)
parser.add_argument('-o', '--output', dest="output",
                    help="Dereplicated fasta filename. Default = 'per_sample_derep.fasta", default='per_sample_derep.fasta')
args = parser.parse_args()

def main():
    print("\n--------------------\n")
    print("Running virome_per_sample_derep.py\n")
    ### Import contig_IDs from summary table
    contig_IDs = list(set(pd.read_csv(args.summary_table, sep='\t')['contig_ID'].tolist()))
    # extract VIBRANT prophage list (contig_IDs containing 'fragment')
    vibrant_prophage = [i for i in contig_IDs if 'fragment' in i]
    # extract VirSorter2 prophage list (contig_IDs containing 'partial')
    vsort2_prophage = [i for i in contig_IDs if 'partial' in i]
    # Remove any retained full contigs from contig_IDs if excised versions have been included by VirSorter2 or VIBRANT (i.e. only keep excised/trimmed versions)
    if args.vibrant_fasta is not None:
        for item in vibrant_prophage:
            ID_trim = re.sub(r'_fragment.*', r'', item)
            # filter out any full contigs of prophage (keeping 'fragment' version)
            contig_IDs = [ x for x in contig_IDs if x != ID_trim ]
    if args.vsort2_fasta is not None:
        for item in vsort2_prophage:
            ID_trim = re.sub(r'\|\|\d+_partial.*', r'', item)
            # filter out any full contigs of prophage (keeping 'partial' version)
            contig_IDs = [ x for x in contig_IDs if x != ID_trim ]
    ### Extract subset of putative viral contigs from original assembly fasta file.
    # Read assmebly fasta, extract contigs found in 'contig_IDs' set, and write to new per_sample_derep.fasta
    with open(args.assembly_fasta, 'r') as read_fasta:
        with open(args.output, 'w') as write_fasta:
            for name, seq in SimpleFastaParser(read_fasta):
                if name in contig_IDs:
                    write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")
    # If VIBRANT fasta file provided, append any prophage identified by vibrant   
    if args.vibrant_fasta is not None:
        with open(args.vibrant_fasta, 'r') as read_fasta:
            with open(args.output, 'a') as write_fasta:
                for name, seq in SimpleFastaParser(read_fasta):
                    if name in vibrant_prophage:
                        write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")
    # If VirSorter2 fasta file provided, append any prophage identified by VirSorter2   
    if args.vsort2_fasta is not None:
        with open(args.vsort2_fasta, 'r') as read_fasta:
            with open(args.output, 'a') as write_fasta:
                for name, seq in SimpleFastaParser(read_fasta):
                    if name in vsort2_prophage:
                        write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")
    # END
    print("Output:\n")
    print(args.output + " : Dereplicated fasta file\n")
    print("Completed virome_per_sample_derep.py\n")
    print("--------------------\n")


if __name__ == '__main__':
    main()

