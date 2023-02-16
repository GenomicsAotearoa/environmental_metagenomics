# Introduction

The examples below are based on sequencing of 12 isolates, multiplexed during ONT library prep (each tagged with a unique barcode). Long-read sequencing data was generated on an Oxford Nanopore GridION, with basecalling based on high accuracy basecalling (HAC; Q9). In this example, basecalling was done in real time during sequencing. However, you can also opt to perform basecalling separately on the raw data after the fact. In our tests, comparable (generally identical) steps were also appropriate for processing data based on super accuracy basecalling model (SUP; Q10).

> **HAC vs SUP**<br>
> Using HAC data can increase coverage (compared with SUP data), but at the expense of a lower quality threshold for individual reads. Coverage depth can be an important factor in the completeness and accuracy of genome assemblies downstream, and so can be important to consider when choosing between HAC and SUP data.

# Data preparation

### Concatenate read chunks per barcode

Outputs from minION/gridION basecallers come in multiple chunks for each barcode. These chunks need to be concatenated to create one fastq file per barcode. Let's assume a directory structure `fastq_files/barcode*/*.fastq.gz`. You can concatenate chunks per barcode:

```bash
# Output directories
mkdir -p concat_barcodes

for barcode_path in fastq_files/barcode*; do
  barcode=$(basename ${barcode_path})
  cat ${barcode_path}/*.fastq.gz > concat_barcodes/${barcode}.fastq.gz
done  
```

> **Reads assigned to unused barcodes**<br>
> Unused barcodes may have a few incorrectly assigned reads due to cross-talk, resulting in a low error rate. It's important to keep in mind that there could be some inaccurate reads in your sample data sets, but these can be identified downstream based on differences in coverage between the correct sequences and those from cross-talk. For now, you can disregard reads in the unused barcode directories.

### Concatenation strategies depending on experimental design

Based on your experimental design, you may wish to further concatenate sequences based on some  criteria.

1. Replicate samples sequenced across multiple barcodes and/or runs to increase data volume
2. Replicate isolates sequenced across multiple barcodes in the same run

So long as you know which barcodes correspond to which sample/isolate, you can further pool/concatenate the relevant fastq files prior to sequence QC and assembly. For example, lets assume you have a file that maps samples/isolates to their barcodes called `sample2barcode.csv`:

```txt
sampleA,barcode001
sampleA,barcode002
sampleA,barcode003
sampleB,barcode004
sampleB,barcode005
sampleC,barcode006
sampleC,barcode007
sampleC,barcode008
...
...
```

You could concatenate files like so:

```bash
# Output directory
mkdir -p concat_sample

while read -r line; do
  sample=$(echo $line | cut -d ',' -f 1)
  barcode=$(echo $line | cut -d ',' -f 2)
  cat ${barcode}.fastq.gz >> ${sample}.fastq.gz
done < sample2barcode.csv
```

Alternatively, you can process sequences belonging to each barcode independently and generate multiple assemblies from the same sample/isolate then dereplicate them later.

<br>

# Quality control of reads

There are several options to assess read quality statistics. In this example we will generate basic read quality stats via NanoStat. 

Note: 
- It may take a few mins per sample to run this. If connection drop outs are an issue, you could run this remotely via, e.g., slurm or tmux
- If you get very fragmented assemblies downstream, you opt to apply a filter here (e.g. length or quality filter) to remove some of the poorer reads, which may improve assemblies downstream. Some options for this include nanofilt or filtlong.

### Run Nanostat on data for all 12 isolates

```bash
# Working directory
cd /working/dir
mkdir -p 1.ONT_data_HAC_concatenated/1.concat_replicates/NanoStat

# Load modules
module purge
module load NanoStat/1.5.0-gimkl-2020a-Python-3.8.2

# Run NanoStat on each barcode dataset
for i in {1..12}; do
    NanoStat -t 8 --tsv \
    -n 1.ONT_data_HAC_concatenated/1.concat_replicates/NanoStat/isolate_${i}_NanoStat.tsv \
    --fastq 1.ONT_data_HAC_concatenated/1.concat_replicates/isolate_${i}.fastq.gz
done
```

### Merge Nano stat results into one summary table[^1]

```bash
# Working directory
cd /working/dir

# Load python
module purge
module load Python/3.8.2-gimkl-2020a
python3

### Import required libraries
import pandas as pd
import numpy as np

# Compile NanoStat results
results_list = []
for i in range (1,13):
    tmp_df = pd.read_csv('1.ONT_data_HAC_concatenated/1.concat_replicates/NanoStat/isolate_'+str(i)+'_NanoStat.tsv', sep='\t')
    tmp_df.index=['isolate_'+str(i)]*len(tmp_df)
    tmp_df = tmp_df.pivot(columns='Metrics', values='dataset')
    results_list.append(tmp_df)

# Generate summary table and write out
results_df = pd.concat(results_list, axis=0)
results_df.index.name = 'isolateID'
results_df.to_csv('1.ONT_data_HAC_concatenated/1.concat_replicates/NanoStat/summary_table_NanoStat.tsv', sep='\t')

quit()
```

<br>

# Assembly

There a several options for assemblers for long read data. In this example, we will use Flye.

Useful outputs:
- `assembly.fasta`: draft assembly
- `tail flye.log`: basic assembly statistics
- `assembly_info.txt`: more information on each contig of the assembly

Note:
- Flye has an option to account for uneven coverage (e.g. metagenomics data or data generated from non-pure isolates): --meta. This can optionally be included in case there are more than one organism in what might be presumed to be a pure isolate (or if the originating sample was already known to be mixed (i.e. more than one genome in the sample))
- Flye has three nano options
  - `nano_corr`: Assumes data are pre-cleaned via, e.g., illumina reads, prior to assembly
  - `nano_raw`: Assumes older poor quality reads, but can also use on HAC data
  - `nano_hq`: Good for SUP-called reads, but can also be used for HAC reads
    - When working with HAC data, this option might return more split contigs. However, preliminary testing in our group did not suggested an appreciable difference between nano_raw and nano_hq with HAC data. You may wish to try both with your data and assess the assemblies generated.

### Run Flye

Run Flye assemblies for each barcode via slurm array (#SBATCH --array=1-12)

Note:
- Running with --meta flag and --nano-hq setting

```bash


#!/bin/bash
#SBATCH -A your_project_account
#SBATCH -J 2_assembly
#SBATCH --time 00:45:00
#SBATCH --array=1-12
#SBATCH --ntasks=1
#SBATCH --mem 10GB
#SBATCH --cpus-per-task=16
#SBATCH -e 2_assembly_%a.err
#SBATCH -o 2_assembly_%a.out

# Working directory
cd /working/dir
mkdir -p 2.assembly/1.assembly.flye.nano_hq

# load module
module purge
module load Flye/2.9-gimkl-2020a-Python-3.8.2

# Run assemblies for each barcode
srun flye --meta -t 16 \
--nano-hq 1.ONT_data_HAC_concatenated/1.concat_replicates/isolate_${SLURM_ARRAY_TASK_ID}.fastq.gz \
--o 2.assembly/1.assembly.flye.nano_hq/isolate_${SLURM_ARRAY_TASK_ID}.assembly
```

Note:
- Some assemblies may fail (e.g. where no contigs are assembled).
- This may be due to, for example; poor coverage for this isolate; insufficient data (similar to poor coverage); mixed genomes in a single sample that the assembler is struggling to resolve.
- Take a note of failed isolateIDs, as this will be useful to know to exclude these isolates downstream.

### Generate summary table of assembly_info.txt outputs

Merge Flye assembly_info.txt outputs into one summary table[^1]

Note:
- The code below includes a try-execpt clause to allow for cases where assembly failed for some isolates

```bash
# Working directory
cd /working/dir

# Load python
module purge
module load Python/3.8.2-gimkl-2020a
python3

### Import required libraries
import pandas as pd
import numpy as np

# Compile flye assembly_info results
results_list = []
for i in range (1,12):
    try:
        tmp_df = pd.read_csv('2.assembly/1.assembly.flye.nano_hq/isolate_'+str(i)+'.assembly/assembly_info.txt', sep='\t').rename({'#seq_name': 'contigID'}, axis=1)    
    except FileNotFoundError:
        continue
    tmp_df.index = ['isolate_'+str(i)]*len(tmp_df)
    tmp_df = tmp_df.rename_axis('isolateID').reset_index()
    results_list.append(tmp_df)

# Generate summary table and write out
results_df = pd.concat(results_list, axis=0)
results_df.to_csv('2.assembly/1.assembly.flye.nano_hq/summary_table_flye_nano_hq_assembly_info.tsv', sep='\t', index=False)

quit()
```

<br>

# Long read polishing

Long-read polishing takes the assemblies that were generated above and maps reads to them and checks and revises them. There are several options for this step, such as racon and medaka. In the example below, we will use medaka.

medaka is Oxford Nanopore's in-house assembly polisher. It takes sequences, assembly files, and the model used for basecalling, and generates polished assembly files.

Note:
- to find model options, run: medaka tools list_models, and then select the applicable model for your basecalling on this specific run.

### Run medaka polishing of Flye assembly

Run medaka polishing and then copy all polished assembly files into one output directory

NOTE:
- In this example, the applicable model was r941_min_hac_g507. Update -m r941_min_hac_g507 if the basecalling model was different for your run.
- If some samples failed the assembly step, you can exclude these from the slurm array in this script, as they will only result in failed runs here.
  - e.g. The example below assumes that samples 2 and 7 failed to assemble any contigs.

```bash
#!/bin/bash
#SBATCH -A your_project_account
#SBATCH -J 2_assembly_polish
#SBATCH --time 01:00:00
#SBATCH --array=1,3-6,8-12
#SBATCH --ntasks=1
#SBATCH --mem 10GB
#SBATCH --cpus-per-task=8
#SBATCH -e 2_assembly_polish_%a.err
#SBATCH -o 2_assembly_polish_%a.out

# Working directory
cd /working/dir
mkdir -p 2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_files

# load module
module purge
module load medaka/1.6.0-Miniconda3-4.12.0

# Run medaka
medaka_consensus \
-m r941_min_hac_g507 \
-i 1.ONT_data_HAC_concatenated/1.concat_replicates/isolate_${SLURM_ARRAY_TASK_ID}.fastq.gz \
-d 2.assembly/1.assembly.flye.nano_hq/isolate_${SLURM_ARRAY_TASK_ID}.assembly/assembly.fasta \
-o 2.assembly/2.assembly.flye.nano_hq.LR_polished/isolate_${SLURM_ARRAY_TASK_ID}.assembly.polished

# Copy assembly file to final output directory
cp 2.assembly/2.assembly.flye.nano_hq.LR_polished/isolate_${SLURM_ARRAY_TASK_ID}.assembly.polished/consensus.fasta \
2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_files/isolate_${SLURM_ARRAY_TASK_ID}.consensus.fasta
```

### Evaluating polished assemblies

It always pays to assess the quality of final assemblies. There as numerous tools to do this and metrics you can assess. Here we will get basic stats via contig counts, contig lengths, and NL50 stats via stats.sh from the BBMap suite of tools. We can then use a bash script to generate summary table from BBMap's stats.sh output.

NOTE:
- The script below runs in a loop over all 12 isolateID numbers, but some isolates may not have files available here if they failed to assemble any contigs in the assembly step (in this example, isolate_2 and isolate_7). We can ignore the following error message: Exception in thread "main" java.lang.RuntimeException: Input file does not appear to be valid (assuming they match the samples that failed to assemble any contigs in the assembly step above)

```bash
# Working directory
cd /working/dir
mkdir -p 2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_stats

# load module
module purge
module load BBMap/38.95-gimkl-2020a

# Run stats on all isolate assemblies
for i in {1..12}; do
    stats.sh in=2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_files/isolate_${i}.consensus.fasta \
    > 2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_stats/isolate_${i}.stats
done

# Extract key stats from outputs and compile into single summary table
echo -e "isolateID\tcontig_count\tcontig_length_max\tcontig_length_total\tcontig_NL50" \
> 2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_stats/summary_table_LRpolished_assembly_stats.tsv
for i in {1..12}; do
    infile="1.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_stats/isolate_${i}.stats"
    isolateID="isolate_${i}"
    contig_count=$(sed -n -e 's/Main genome contig total:[[:space:]]\+//p' ${infile})
    contig_length_max=$(sed -n -e 's/Max contig length:[[:space:]]\+//p' ${infile})
    contig_length_total=$(sed -n -e 's/Main genome contig sequence total:[[:space:]]\+\(.*MB\).*/\1/p' ${infile})
    contig_NL50=$(sed -n -e 's/Main genome contig N\/L50:[[:space:]]\+//p' ${infile})
    echo -e "${isolateID}\t${contig_count}\t${contig_length_max}\t${contig_length_total}\t${contig_NL50}" \
    >> 2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_stats/summary_table_LRpolished_assembly_stats.tsv
done
```

# Short read polishing

Historically, Nanopore long-read sequencing has suffered from high error rates, and methods were developed to polish long-read assemblies via adding high quality short read (such as Illumina HiSeq) data as well. Some assembly algorithms have also been developed that take both long read and short read data together when generating assemblies.

Recent advances have considerably increased Nanopore basecalling quality, and ensuring a high sequencing coverage across the genome also likely helps mitigate some of the remaining errors.

In one small pilot test within our group, checkM completeness and contamination stats were comparable between only long-read polished isolate data and both long- and short-read polished data. Although, the quality and consistency at the individual base level was not tested, and how important this is for your study may depend on what you're asking of the data downstream.

We will not cover the process of short-read polishing in these docs, but it is worth considering for your data and study question whether you wish to include this step here. Ultimately, if you have high quality short-read data available, or it is easily obtainable, it is likely to be beneficial to incorporate them here.

# Assembly post-processing

### Add isolateIDs to contig headers

It can be useful for downstream applications (particularly where assembled contigs from all sequenced isolates are analysed together) to add isolateIDs into the headers of contigIDs for each isolate assembly.

The script below uses sed to do this in place on the polished_assembly_files/. (If you wish to keep these original files unedited, you can remove the -i flag and instead read the result into a new file via >)

```bash
cd /working/dir

# nano-hq assemblies
for i in {1..12}; do
    for file in 2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_files/isolate_${i}.consensus.fasta; do
        sed -i -e "s/>/>isolate_${i}_/g" ${file}
    done
done
```

### Generate assembly2contig_lookupTable

A lookup table mapping all contig IDs to assembly IDs for all assemblies can be useful in some cases for downstream applications.

If this is of use, you can use the python code below to generate assembly2contig_lookupTable.txt

```bash
# Working directory
cd /working/dir

# Load python
module purge
module load Python/3.8.2-gimkl-2020a
python3

### Import required libraries
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
import os

# Compile genome2contig lookup tables
with open('2.assembly/2.assembly.flye.nano_hq.LR_polished/assembly2contig_lookupTable.tsv', 'w') as write_file:
    # header
    write_file.write('assemblyID\tcontigID\n')
    # Prok bins
    assemblyfiles_directory = os.fsencode('2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_files')
    for file in os.listdir(assemblyfiles_directory):
        filename = os.fsdecode(file)
        assemblyID = os.path.splitext(filename)[0]
        with open('2.assembly/2.assembly.flye.nano_hq.LR_polished/polished_assembly_files/' + filename, 'r') as read_fasta:
            for name, seq in SimpleFastaParser(read_fasta):
                write_file.write(assemblyID + '\t' + name + '\n')

quit()
```

[^1]: This will be written into a script on a later date
