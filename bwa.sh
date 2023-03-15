#!/usr/bin/env bash

fasta=$1
r1=$2
r2=$3
unpaired=$4
outdir=$5
threads=$6

mkdir -p $outdir

#bwa index $fasta

minimap2 -ax sr -F 800 -t $threads $fasta $unpaired > $outdir"/aln-se.sam"

minimap2 -ax sr -F 800 -t $threads $fasta $r1 $r2 > $outdir"/aln-pe.sam"

#bwa mem -t $threads $fasta $r1 $r2 > $outdir"/aln-pe.sam"
#bwa mem -t $threads $fasta $unpaired > $outdir"/aln-se.sam"

sam_to_bam.sh $outdir"/aln-pe.sam" $threads
sam_to_bam.sh $outdir"/aln-se.sam" $threads

pilon \
--genome $fasta \
--frags $outdir"/aln-pe.bam" \
--unpaired $outdir"/aln-se.bam" \
--outdir $outdir"/pilon" \
--threads $threads \
--changes \
--vcf
