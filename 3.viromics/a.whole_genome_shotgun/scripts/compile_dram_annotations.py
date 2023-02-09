#!/usr/bin/env python

'''compile_dram_annotations.py

Compile annotations.tsv, trnas.tsv, and rrnas.tsv output files from DRAM subset runs (e.g. from manually parallelised runs)


Required parameters:

--input (-i) : Path to directory containing output directories for each DRAM(-v) subset

Optional parameters:

--output (-o) : Optional prefix for output filenames. Can include directory path (e.g. 'path/to/output/collated_dram').

'''

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import os.path
from glob import glob

parser = ArgumentParser()
parser.add_argument('-i', '--input', dest="inpath",
                    help="Path to directory containing output directories for each DRAM(-v) subset", required=True)
parser.add_argument('-o', '--output', dest="outprefix",
                    help="Optional prefix for output filenames. Can include directory path (e.g. 'path/to/output/collated_dram').", default='')
args = parser.parse_args()

def main(): 
    print("\n--------------------\n")
    print("Running compile_dram_annotations.py\n")
    # Loop through the three file types (annotations, rrnas, trnas)
    for file in ['annotations', 'rrnas', 'trnas']:   
        # List of dram subset output directories
        dataframe_paths = sorted(glob(args.inpath+'/*/'+file+'.tsv'))
        # Concatenate subsets (if no objects to concatenated, output empty file)
        try:
            merged_subsets = pd.concat([pd.read_csv(i, sep='\t', index_col=0) for i in dataframe_paths])
        except ValueError:
            merged_subsets = pd.DataFrame()
        # Write out concatenated dataframe
        merged_subsets.to_csv(args.outprefix+file+'.tsv', sep='\t', index=True)
    print("Completed compile_dram_annotations.py\n")
    print("--------------------\n")


if __name__ == '__main__':
    main()

