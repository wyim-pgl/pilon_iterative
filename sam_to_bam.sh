#!/usr/bin/env bash


sam=$1
threads=$2
bam=${sam/.sam/.bam}

cat $sam | samtools view -@ $threads -h  -f 2 -F 512 | samtools sort -@ $threads -m 10G -l0 > $bam

samtools index -@ $threads $bam
