# Viral metagenomics

This workflow was developed within the group of Associate Professor Kim Handley (University of Auckland) under the Environmental Metagenomics branch of [Genomics Aotearoa](https://www.genomics-aotearoa.org.nz/) ([GitHub](https://github.com/GenomicsAotearoa/)).

## Data origin and assumptions

Data used to inform the preparation of this workflow was generated from multiple water samples sequenced on the Illumina HiSeq platform for the study of DNA viruses with a focus on phages. The process is comparable for data generated from Oxford Nanopore long read sequencing of isolates, but take note of the following:

- DNA samples in the study used to generate the [long read sequence processing workflow](../../1.sequence_processing/4.isolates_longread.md) were obtained from bacterial isolates grown in liquid media then spun down into a pellet before extraction. For viral identification, note that this process will predominantly target prophages integrated into the host genome at the time of DNA extraction. Lytic viruses are less likely to be caught here. However, in our experience, there may also be some cases where assembled contigs are fully circular viral genomes, which may represent intracellular viruses that are not integrated into the genome (e.g. replicating at the time; or viruses that are intracellular but do not integrate into the genome (akin to a plasmid)), or extracellular viral particles caught within the pellet during centrifugation and DNA extraction.
- It may be preferable to skip the final dereplication step (via `Cluster_genomes.pl`). For individual isolate genomes, you will likely want to retain closely related viruses as distinct genomes rather than clustering them together (especially if they originated from different isolates, or are from the same isolate but one sequence is excised from an integrated prophage, and another identical sequence is a circular genome as this would suggest both integration and replication in the host).

<!--
If you only have a single co-assembly, then modify appropriately to only run on that one data set.

If you have both individual assemblies and a co-assembly (and/or mini-co-assemblies) then modify the scripts appropriately to run on all data sets (you can treat the co-assembly as simply an additional assembly data set; so in this case, if you had nine sample assemblies and one co-assembly, you can run these as if you had 10 individual assemblies).
-->

## Software prerequisites
