#!/bin/bash

# Example:
# delete_fasta_segment input.fasta output.fasta       # deletes 10 bp at pos 1000
# delete_fasta_segment in.fa out.fa 500 5             # deletes 5 bp at pos 500
delete_fasta_segment() {
  local infile="$1"
  local outfile="$2"
  local pos="${3:-1000}"
  local len="${4:-10}"

  if [[ -z "$infile" || -z "$outfile" ]]; then
    echo "Usage: delete_fasta_segment <input.fa> <output.fa> [pos] [len]" >&2
    return 1
  fi

  awk -v pos="$pos" -v len="$len" '
    /^>/ {
      if (seq) {
        pre = substr(seq,1,pos-1)
        suf = substr(seq,pos+len)
        newseq = pre suf
        for(i=1; i<=length(newseq); i+=60)
          print substr(newseq,i,60)
      }
      print                # header
      seq = ""             # reset sequence
      next
    }
    { seq = seq $0 }       # accumulate
    END {
      if (seq) {
        pre = substr(seq,1,pos-1)
        suf = substr(seq,pos+len)
        newseq = pre suf
        for(i=1; i<=length(newseq); i+=60)
          print substr(newseq,i,60)
      }
    }
  ' "$infile" > "$outfile"
}

if true; then
    echo "...creating reference and cancer"
    REF=~/ref/Homo_sapiens_assembly38.fasta
    samtools faidx $REF chr22:23178508-23320037 | sed '1s/^>.*$/>bcr/' > fullbcr.fa
    samtools faidx $REF chr9:130711043-130889675 | sed '1s/^>.*$/>abl/' > fullabl.fa
    samtools faidx $REF chr17:7666421-7689490 | sed '1s/^>.*$/>tp53/' > tp53.fa
    samtools faidx $REF chr8:127733434-127744951 | sed '1s/^>.*$/>myc/' > myc.fa

    ## create an indel cancer genome
    delete_fasta_segment myc.fa myc_del.fa 1000 10

    ## create the reference
    cat fullbcr.fa fullabl.fa tp53.fa myc.fa > tiny.fasta
    rm fullbcr.fa fullabl.fa tp53.fa myc.fa

    ## index the reference
    /usr/bin/samtools faidx tiny.fasta
    bwa index tiny.fasta

    ## create a rearrangement cancer genome
    { 
	printf ">BCRABL\n"
	samtools faidx "$REF" chr22:23220950-23255836 \
	    | tail -n +2
	samtools faidx "$REF" chr9:130855888-130872542 \
	    | tail -n +2
    } > BCRABL.fa

    ## create the final reference genome
    cat BCRABL.fa myc_del.fa > cancer.fa
fi

## simulate normal reads
wgsim -e 0.002 -1 150 -2 150 -d 300 -s 30 -N 100000 -S 1337 tiny.fasta sim1.fq sim2.fq

## simulate tumor reads
wgsim -e 0.002 -1 150 -2 150 -d 300 -s 30 -N 10000 -S 1337 cancer.fa sim1_bcr.fq sim2_bcr.fq

cat sim1.fq sim1_bcr.fq > sim1.fastq
cat sim2.fq sim2_bcr.fq > sim2.fastq

## align
bwa mem -R '@RG\tID:RG1\tSM:sample1\tPL:ILLUMINA' tiny.fasta sim1.fastq sim2.fastq | \
  /usr/bin/samtools view -bS - | \
  /usr/bin/samtools sort -o sim.sorted.bam -
/usr/bin/samtools index sim.sorted.bam

## split is at bcr:77330
## del is at myc:1000
