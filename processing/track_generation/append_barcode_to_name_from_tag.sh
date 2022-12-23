#!/bin/bash

input_BAM=$1

samtools view -H $input_BAM > ${input_BAM%.bam}.barcoded.sam

samtools view $input_BAM | awk '
{
    for (i=1; i<=NF; i++){
        if ($i~/^UB/){
            a=$i
            }
        }
    read_name=$1; $1=""; print read_name "_" substr(a,6,21),$0
}' >> ${input_BAM%.bam}.barcoded.sam



