# Workflow: GC, AA, and pI

Workflow for compiling GC content, amino acid proportions, and isoelectric points of predicted genes and proteins, as used in the study *'Genomic markers show the proliferation of DNA viruses is constrained by adaptation to environmental niche'*.

Method:

- Predict genes and protein sequence via DRAMv (prodigal-gv)
- Calculate animo acid proportions per predicted protein
- Generate summaries for AA types including acidic, basic, polar, nonpolar, and charged amino acids (defined as per pepstats (EMBOSS v6.6.0))
- Identify gene GC content
- Output summary table of results
- Also calculate protein isolectric points (pI)

# Software dependencies

- This workflow was developed and tested in Python v3.8.2, and requires the following python libraries: pandas, numpy, re, os, Bio
