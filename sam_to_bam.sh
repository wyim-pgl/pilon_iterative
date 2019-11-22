#!/usr/bin/env bash


sam=$1
threads=$2
bam=${sam/.sam/.bam}

cat $sam | samtools view -@ $threads -h -F4 | samtools sort -@ $threads -m 10G -l0 > $bam

samtools index -@ $threads $bam