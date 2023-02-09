#!/usr/bin/env python

'''summarise_viral_contigs.py

Outputs a summary table of results from tools used to identify viral contigs in WGS or WTS data. 

- For WTS data, allows for inclusion of results from: blastx rdrp; blastx nr; blastn nt; vibrant; virsorter2; deepvirfinder
- For WGS data, allows for inclusion of results from: vibrant; virsorter2; deepvirfinder


Optional parameters: 

--blastn_nt (-n) : Output from blastn search against NCBI nt database (NOTE: the following output format is expected: 6 qseqid sacc salltitles pident length evalue sskingdoms).
--blastx_nr (-x) : Output from diamond blastx search against NCBI nr database (NOTE: the following output format is expected: 6 qseqid qlen sseqid stitle pident length evalue sskingdoms).
--blastx_rdrp (-r) : Output from diamond blastx search against custom viral rdrp database (NOTE: the following output format is expected: qseqid qlen sseqid stitle pident length evalue).
--vibrant (-v) : Summary results file output from vibrant: VIBRANT_summary_results_SampleX.tsv.
--virsorter2 (-s) : Summary results file output from VirSorter2: SampleX-final-viral-score.tsv.
--deepvirfinder (-d) : Summary results file output from deepvirfinder (filtered by score, p.adj value, and sequences matching Eukaryota (via Kraken2)): SampleX.dvfpred_filtered.txt
--out_prefix (-o) : Output summary table prefix (suffixes = .txt and _VIRUSES.txt). Default = 'viral_contigs_summary_table'.
'''

from argparse import ArgumentParser
import os
import pandas as pd
import numpy as np
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = ArgumentParser()
parser.add_argument('-n', '--blastn_nt', dest="nt_infile",
                    help="Output from blastn search against NCBI nt database (NOTE: the following output format is expected: 6 qseqid sacc salltitles pident length evalue sskingdoms)", 
                    default=None)
parser.add_argument('-x', '--blastx_nr', dest="nr_infile",
                    help="Output from diamond blastx search against NCBI nr database (NOTE: the following output format is expected: 6 qseqid qlen sseqid stitle pident length evalue sskingdoms)",
                   default=None)
parser.add_argument('-r', '--blastx_rdrp', dest="rdrp_infile",
                    help="Output from diamond blastx search against custom viral rdrp database (NOTE: the following output format is expected: qseqid qlen sseqid stitle pident length evalue)",
                   default=None)
parser.add_argument('-v', '--vibrant', dest="vibrant_infile",
                    help="Summary results file output from vibrant: VIBRANT_summary_results_SampleX.tsv",
                   default=None)
parser.add_argument('-s', '--virsorter2', dest="vsort2_infile",
                    help="Summary results file output from VirSorter2: SampleX-final-viral-score.tsv",
                   default=None)
parser.add_argument('-d', '--deepvirfinder', dest="dvf_infile",
                    help="Summary results file output from deepvirfinder (filtered to remove Eukarotic contigs): SampleX.dvfpred_filtered.txt",
                   default=None)
parser.add_argument('-o', '--out_prefix', dest="out_prefix",
                    help="Output summary table prefix (suffixes = .txt and _VIRUSES.txt). Default = viral_contigs_summary_table", default='viral_contigs_summary_table')
args = parser.parse_args()

def main():
    print("\n--------------------\n")
    print("Running summarise_viral_contigs.py\n")
    # Set up emtpy summary_table DataFrame
    summary_table = pd.DataFrame(data={'contig_ID': []})
    # Add NT blastn results to summary_table (if provided)
    if args.nt_infile is not None:
        nt_in = pd.read_csv(args.nt_infile, 
                            sep='\t', 
                            names=["nt_contig_ID", "nt_sacc", "nt_stitle", "nt_pident", "nt_hit_length", "nt_evalue", "nt_Kingdom"])
        nt_in['contig_ID'] = nt_in['nt_contig_ID']
        summary_table = pd.merge(summary_table, nt_in, left_on="contig_ID", right_on="contig_ID", how='outer')
    # Add NR diamond blastx results to summary_table (if provided)
    if args.nr_infile is not None:
        nr_in = pd.read_csv(args.nr_infile, 
                            sep='\t', 
                            names=["nr_contig_ID", "length", "nr_sseqid", "nr_stitle", "nr_pident", "nr_hit_length", "nr_evalue", "nr_Kingdom"])
        nr_in.drop(columns='nr_sseqid', inplace=True)
        nr_in['contig_ID'] = nr_in['nr_contig_ID']
        summary_table = pd.merge(summary_table, nr_in, left_on="contig_ID", right_on="contig_ID", how='outer', suffixes=('', '_DROP')).filter(regex='^(?!.*_DROP)')
    # Add RdRP diamond blastx results to summary_table (if provided)
    if args.rdrp_infile is not None:
        rdrp_in = pd.read_csv(args.rdrp_infile, 
                              sep='\t', 
                              names=["rdrp_contig_ID", "length", "rdrp_sseqid", "rdrp_stitle", "rdrp_pident", "rdrp_hit_length", "rdrp_evalue"])
        rdrp_in.drop(columns='rdrp_sseqid', inplace=True)
        rdrp_in['contig_ID'] = rdrp_in['rdrp_contig_ID']
        summary_table = pd.merge(summary_table, rdrp_in, left_on="contig_ID", right_on="contig_ID", how='outer', suffixes=('', '_DROP')).filter(regex='^(?!.*_DROP)')
    # Add VIBRANT results to summary_table (if provided)
    if args.vibrant_infile is not None:
        vibrant_in = pd.read_csv(args.vibrant_infile, sep='\t')
        vibrant_in.drop(columns=['KEGG int-rep', 'KEGG zero', 'Pfam int-rep', 'Pfam zero', 'VOG redoxin', 'VOG rec-tran', 'VOG int', 'VOG RnR', 'VOG DNA', 'KEGG restriction check', 'KEGG toxin check', 'VOG special', 'annotation check', 'p_v check', 'p_k check', 'k_v check', 'k check', 'p check', 'v check', 'h check'], inplace=True)
        vibrant_in.columns = ['vibrant_contig_ID', 'vibrant_total_genes', 'vibrant_KEGG_genes', 'vibrant_KEGG_v_score', 'vibrant_Pfam_genes', 'vibrant_Pfam_v_score', 'vibrant_VOG_genes', 'vibrant_VOG_v_score']
        vibrant_in['contig_ID'] = vibrant_in['vibrant_contig_ID']
        summary_table = pd.merge(summary_table, vibrant_in, left_on="contig_ID", right_on="contig_ID", how='outer')
    # Add VirSorter2 results to summary_table (if provided)
    if args.vsort2_infile is not None:
        vsort2_in = pd.read_csv(args.vsort2_infile, sep='\t')
        vsort2_in.drop(columns=['dsDNAphage', 'NCLDV', 'RNA', 'ssDNA', 'lavidaviridae', 'length'], inplace=True)
        vsort2_in.columns = ['vsort2_contig_ID', 'vsort2_max_score', 'vsort2_max_score_group', 'vsort2_hallmark_genes', 'vsort2_viral_component', 'vsort2_cellular_component']
        vsort2_in['contig_ID'] = vsort2_in['vsort2_contig_ID'].replace(to_replace=r'\|\|full', value=r"", regex=True)
        vsort2_in['contig_ID'] = vsort2_in['contig_ID'].replace(to_replace=r'\|\|lt2gene', value=r"", regex=True)
        summary_table = pd.merge(summary_table, vsort2_in, left_on="contig_ID", right_on="contig_ID", how='outer')
    # Add DeepVirFinder results to summary_table (if provided)
    if args.dvf_infile is not None:    
        dvf_in = pd.read_csv(args.dvf_infile, sep='\t')
        dvf_in.drop(columns=['len', 'pvalue', 'qvalue'], inplace=True)
        dvf_in.columns = ['dvfind_contig_ID', 'dvfind_score', 'dvfind_p_adj']
        dvf_in['contig_ID'] = dvf_in['dvfind_contig_ID']
        summary_table = pd.merge(summary_table, dvf_in, left_on="contig_ID", right_on="contig_ID", how='outer')
    # Shift 'contig_ID' to start of summary_table
    summary_table = summary_table[ ['contig_ID'] + [ col for col in summary_table.columns if col != 'contig_ID' ] ]
    # Output summary table
    # If blast search outputs provided, output subset summary table, selecting only matches relevant to viruses
    ### NOTE: currently hard-coded assuming all working options provided (blast searches x3, vibrant, vsort2 // excl. dvfind)
    if args.rdrp_infile is not None and args.nt_infile is not None and args.nr_infile is not None and args.vibrant_infile is not None and args.vsort2_infile is not None:
        subset_table = summary_table[(summary_table['nt_stitle'].str.contains("virus|virid|viral").fillna(False)) | (summary_table['nt_Kingdom'].str.contains("Virus").fillna(False)) | (summary_table['nr_stitle'].str.contains("virus|virid|viral").fillna(False)) | (summary_table['nr_Kingdom'].str.contains("Virus").fillna(False)) | (summary_table['rdrp_contig_ID'].notna()) | (summary_table['vibrant_contig_ID'].notna()) | (summary_table['vsort2_contig_ID'].notna())]
        subset_table.to_csv(args.out_prefix+'_VIRUSES.txt', sep='\t', index=False)
        summary_table.to_csv(args.out_prefix+'.txt', sep='\t', index=False)
        print("Output:\n")
        print(args.out_prefix+".txt : Output summary table (Full table)")
        print(args.out_prefix+"_VIRUSES.txt : Output summary table (Filtered for viruses)\n")
    else:
        summary_table.to_csv(args.out_prefix+'_VIRUSES.txt', sep='\t', index=False)
        print("Output:\n")
        print(args.out_prefix+"_VIRUSES.txt : Output summary table\n")
    print("Completed summarise_viral_contigs.py\n")
    print("--------------------\n")


if __name__ == '__main__':
    main()


