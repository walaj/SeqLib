#!/bin/bash

#REF=~/Dropbox_HMS/ref/hg38/Homo_sapiens_assembly38.fasta
#samtools faidx $REF chr22:23178508-23320037 | sed '1s/^>.*$/>bcr/' > fullbcr.fa
#samtools faidx $REF chr9:130711043-130889675 | sed '1s/^>.*$/>abl/' > fullabl.fa
#samtools faidx $REF chr17:7666421-7689490 | sed '1s/^>.*$/>tp53/' > tp53.fa
#samtools faidx $REF chr8:127733434-127744951 | sed '1s/^>.*$/>myc/' > myc.fa
#cat fullbcr.fa fullabl.fa tp53.fa myc.fa > tiny.fa

#rm fullbcr.fa fullabl.fa tp53.fa myc.fa
#samtools faidx tiny.fa
#bwa index tiny.fa

# { 
#   printf ">BCRABL\n"
#   samtools faidx "$REF" chr22:23220950-23255836 \
#     | tail -n +2
#   samtools faidx "$REF" chr9:130855888-130872542 \
#     | tail -n +2
# } > BCRABL.fa

## simulate normal reads
wgsim -e 0.002 -1 150 -2 150 -d 300 -s 30 -N 100000 tiny.fa sim1.fq sim2.fq

## simulate tumor reads
wgsim -e 0.002 -1 150 -2 150 -d 300 -s 30 -N 10000 BCRABL.fa sim1_bcr.fq sim2_bcr.fq

cat sim1.fq sim1_bcr.fq > sim1.fastq
cat sim2.fq sim2_bcr.fq > sim2.fastq

## align
bwa mem -R '@RG\tID:RG1\tSM:sample1\tPL:ILLUMINA' tiny.fa sim1.fastq sim2.fastq | \
  samtools view -bS - | \
  samtools sort -o sim.sorted.bam -
samtools index sim.sorted.bam

## split is at bcr:77212
