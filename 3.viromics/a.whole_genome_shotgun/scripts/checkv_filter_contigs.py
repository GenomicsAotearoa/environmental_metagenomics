#!/usr/bin/env python

'''checkv_filter_contigs.py

Filter fasta based on CheckV results: retain only contigs where: (viral gene > 0 OR (viral gene == 0 AND host gene == 0))

- Writes out new filtered fna file and filtered checkv quality summary tsv file


Required parameters:

--checkv_dir_input (-i) : Directory containing checkv output files viruses.fna, proviruses.fna, and quality_summary.tsv

Optional parameters:

--output_prefix (-o) : Prefix for output files. Default = checkv_contigs
'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument("-i", "--checkv_dir_input", dest="checkv_dir", 
                    help="Directory containing checkv output files viruses.fna, proviruses.fna, and quality_summary.tsv",
                    metavar='CheckV output directory', required=True)
parser.add_argument('-o', '--output_prefix', dest="output",
                    help="Prefix for output files. Default = checkv_contigs", default='checkv_contigs')
args = parser.parse_args()

def import_quality_summary():
    # Import checkv quality file
    qual_summary = pd.read_csv(args.checkv_dir+'/quality_summary.tsv', sep='\t')
    return qual_summary

def import_fna_files():
    # Read fasta files from checkv output (viruses.fna and proviruses.fna) into fasta dictionary
    ## In the process, modify provirus contig headers for easier downstream use
    ## set up dictionary
    fasta_dict = {}
    ## read viruses.fna into dictionary
    with open(args.checkv_dir+'/viruses.fna', 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            fasta_dict[name] = seq
    ## read proviruses.fna into dictionary
    with open(args.checkv_dir+'/proviruses.fna', 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            name_edit = re.sub(r'\s', '__checkv_excised_start_', name)
            name_edit = re.sub(r'-', '_end_', name_edit)
            name_edit = re.sub(r'\/', '_len_', name_edit)
            fasta_dict[name_edit] = seq
    return fasta_dict

def checkv_filt(qual_summary, fasta_dict): 
    print('Filtering contigs...')
    # Subset contigs by minimum checkv virus and host gene counts thresholds
    contig_ids_filt = list(qual_summary.loc[(qual_summary['viral_genes'] > 0) | ((qual_summary['viral_genes'] == 0) & (qual_summary['host_genes'] == 0))]['contig_id'])
    print("Total contigs = " + str(len(fasta_dict)))
    print("Retained contigs = " + str(len(contig_ids_filt)))
    # Write out filtered fna file
    print("Writing filtered fna file: " + args.output + "_filtered.fna\n")
    with open(args.output+'_filtered.fna', 'w') as write_fasta:
        for name, seq in fasta_dict.items():
            name_lookup = re.sub(r'_[1-9]__.*', '', name)
            if name_lookup in contig_ids_filt:
                write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")
    # Write out filtered quality summary file
    print("Writing filtered checkv quality_summary: " + args.output + "_filtered_quality_summary.tsv\n")
    qual_summary_filt = qual_summary[qual_summary['contig_id'].isin(contig_ids_filt)]
    qual_summary_filt.to_csv(args.output+'_filtered_quality_summary.tsv', sep='\t')

def main():
    print("\n--------------------\n")
    print("Running checkv_filter_contigs.py\n")
    qual_summary = import_quality_summary()
    fasta_dict = import_fna_files()
    # Filter by checkv viral gene and host gene counts
    checkv_filt(qual_summary, fasta_dict)
    print("Completed checkv_filter_contigs.py\n")
    print("--------------------\n")


if __name__ == '__main__':
    main()


