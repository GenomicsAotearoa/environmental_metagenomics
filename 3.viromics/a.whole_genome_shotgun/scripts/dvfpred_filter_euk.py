#!/usr/bin/env python

'''dvfpred_filter_euk.py

Filter Eukaryota sequences out of DeepVirFinder results
'''

from argparse import ArgumentParser
import os
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument("-d", "--deepvirfinder_results", dest="infile_dvf",
                    help="Unfiltered results file from DeepVirFinder: dvfpred_unfiltered/<sample>_gt1000bp_dvfpred.txt",
                    metavar='DeepVirFinder_unfiltered_results', required=True)
parser.add_argument("-f", "--Euk_fasta", dest="infile_fasta",
                    help="Fasta file of Eukaryota sequences identified by Kraken2. e.g. dvfpred_kraken/<sample>.kraken_Euk.fasta",
                    metavar='Eukaryota_fasta', required=True)
parser.add_argument('-o', '--output', dest="output_dvf_filtered",
                    help="Output filename for filtered DeepVirFinder results", required=True)
args = parser.parse_args()

def main():
    print("\n--------------------\r\n")
    print("Running dvfpred_filter_euk.py for "+args.infile_dvf+"\r\n")
    # Import data into python
    DVF = pd.read_csv(args.infile_dvf, sep='\t')
    # Generate set of contig IDs identified as Eukaryota via Kraken2.
    ## Read output from extract_kraken_reads.py (<sample>.kraken_Euk.fasta), extract sequence IDs to then filter out of dvfpred results.
    DVF_eukID = set()
    with open(args.infile_fasta, 'r') as read_fasta:
        for name, seq in SimpleFastaParser(read_fasta):
            DVF_eukID.add(name)
    # Filter Eukaryota sequence IDs (stored as a set in DVF_eukID) out of DVF results DataFrame
    print("Filtering Eukaryota contigs out of DeepVirFinder predictions...")
    DVF = DVF[~DVF['name'].isin(list(DVF_eukID))]
    # Write filtered DVF results to file
    print("Writing filtered results to file...\r\n")
    DVF.to_csv(args.output_dvf_filtered, sep='\t', index=False)
    print("Output:\r\n")
    print(args.output_dvf_filtered+": DeepVirFinder predictions, filtered to remove sequences matching Eukaryota\r\n")
    print("Completed dvfpred_filter_euk.py for "+args.infile_dvf+"\r\n")
    print("--------------------\r\n")

if __name__ == '__main__':
    main()
