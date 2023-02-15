# MAG construction

## Introduction

Using inherent compositional properties, such as oligonucleotide frequencies and coverage profiles, of the assembled contigs/scaffolds, we can group contigs/scaffolds into "bins". Here, a "bin" (or MAG) essentially a population genome. In general, this is a 3-step process:

1. Multi-tool binning
2. Dereplication
3. Manual refinement

How step 2 is performed in this process depends on your assembly strategy, where if you have:

- One (co-)assembly, you only need to dereplicate across tools/software
- Multiple (co-)assemblies, you need to first dereplicate *across tools/software per assembly*, ***then*** dereplicate *across assemblies*. 

Bin evaluation are also performed at steps 2 and 3.

## Required data

- Quality trimmed reads
- Assembled contigs/scaffolds
