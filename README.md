# PWT:  Pilon With iTerations

## Introduction
PWT (pronounced 'pout') is a wrapper for the genome polishing tool [pilon](https://github.com/broadinstitute/pilon).  It uses a pair of Illumina reads in fastq format as well as a file of unpaired reads to improve the quality of a genome assembly.  PWT will run pilon multiple times until pilon doesn't make changes to the assembly.

## Method:
1) align reads to the assembly
2) run pilon
3) classify the changes
4) if there are changes, return to step 1 using the assembly produced by pilon

## Installation

git clone https://github.com/rehrlich/pwt.git

### Install dependencies:
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda pilon
conda install pandas

By default, pilon doesn't use enough memory for many projects.  Run 'which pilon' to find the pilon file, then change default_jvm_mem_opts to values that are appropriate for your system

### Developed using:
bwa 0.7.17-r1188
samtools 1.9
pilon 1.23
pandas 0.25.3

## Usage:
loop_pilon.py [-h] [-m MAX_ITERATIONS] [-t THREADS] --r1 R1 --r2 R2 -u UNPAIRED -f FASTA_PATH -o OUTDIR -s STRAIN_NAME

optional arguments:
  -h, --help            show this help message and exit
  -m MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        Max number of times to run pilon (default: 10)
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 16)

Required arguments:
  --r1 R1               Path to r1 paired reads in fastq format
  --r2 R2               Path to r2 paired reads in fastq format
  -u UNPAIRED, --unpaired UNPAIRED
                        Path to unpaired reads in fastq format
  -f FASTA_PATH, --fasta_path FASTA_PATH
                        Path to reference fasta
  -o OUTDIR, --outdir OUTDIR
                        Directory for output, must not exist
  -s STRAIN_NAME, --strain_name STRAIN_NAME
                        Nickname for the reference fasta

## Format of outdir/changes.csv:
Each row corresponds to one change made by pilon
pilon_num - the iteration number
variant - the type of change made by pilon.  Each change is explained below under 'Variant Types'
size_change - the change in the length of the sequence.  A positive number means that Pilon increased the size of the sequence
source / destination - the source or destination from the pilon output.  Source refers to the reference fasta; destination refers to the fasta produced by pilon
*_seq - the sequence.  A period means the sequence is empty
*_ctg - the name of the contig
*_start - the start position of the variant
*_end - the end position of the variant

## Variant Types:
Insertion - Pilon added one or more consecutive bases to the reference.  Size change is positive
Deletion - Pilon removed one or more consecutive bases from the reference.  Size change is negative
SNP - Pilon changed one base in the reference to a different base
Complex variant - Any change not listed above
