#!/bin/bash

set -x
set -u
set -e

ulimit -n 500000

source ./SLR_pipeline_parameters.sh

outdir=$output_directory


mkdir -p $outdir


###define input fastqs

threeprime_R1=$input_dir_3prime/*R1_001.fastq.gz
threeprime_R2=$input_dir_3prime/*R3_001.fastq.gz
threeprime_I1=$input_dir_3prime/*R2_001.fastq.gz


internal_R1=$input_dir_internal/*R1_001.fastq.gz
internal_R2=$input_dir_internal/*R3_001.fastq.gz
internal_I1=$input_dir_internal/*R2_001.fastq.gz


fiveprime_R1=$input_dir_5prime/*R1_001.fastq.gz
fiveprime_R2=$input_dir_5prime/*R3_001.fastq.gz
fiveprime_I1=$input_dir_5prime/*R2_001.fastq.gz



#export compress_fastq_threads="8"
: '
~/analysis_10/SLR/tools/single_cell_toolkit/barcode_10x_scatac_fastqs.sh \
	 $threeprime_R1 \
         $threeprime_I1 \
         $threeprime_R2 \
         $input_dir_3prime/threeprime_barcoded \
         false \
         false \
         ":" \
	 "igzip"

~/analysis_10/SLR/tools/single_cell_toolkit/barcode_10x_scatac_fastqs.sh \
         $internal_R1 \
         $internal_I1 \
       	 $internal_R2 \
         $input_dir_internal/internal_barcoded \
         false \
         false \
	 ":" \
         "igzip"

~/analysis_10/SLR/tools/single_cell_toolkit/barcode_10x_scatac_fastqs.sh \
         $fiveprime_R1 \
         $fiveprime_I1 \
         $fiveprime_R2 \
         $input_dir_5prime/fiveprime_barcoded \
         false \
         false \
	 ":" \
         "igzip"
'


#Align 3prime reads with 10X barcodes

#threeprime_R1_barcoded=${threeprime_R1%.fastq.gz}.barcoded.fastq.gz
#threeprime_R2_barcoded=${threeprime_R2%.fastq.gz}.barcoded.fastq.gz

threeprime_R1_barcoded=$input_dir_3prime/threeprime_barcoded_R1.fastq.gz
threeprime_R2_barcoded=$input_dir_3prime/threeprime_barcoded_R2.fastq.gz

: '
STAR \
      --genomeDir $reference \
      --genomeLoad LoadAndKeep \
      --limitBAMsortRAM 50000000000 \
      --runThreadN 42 \
      --outSAMattrRGline ID:3prime \
      --soloType CB_UMI_Simple \
      --readFilesIn $threeprime_R2_barcoded $threeprime_R1_barcoded \
      --soloCBstart 1 \
      --soloCBlen 16 \
      --soloUMIstart 17 \
      --soloUMIlen 12 \
      --soloBarcodeReadLength 0 \
      --soloCBwhitelist $BC_whitelist  \
      --readFilesCommand zcat \
      --outSAMattributes GX GN CB CR CY UB UR UY \
      --outSAMtype BAM SortedByCoordinate \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIdedup 1MM_CR \
      --soloFeatures GeneFull Gene SJ \
      --soloMultiMappers Unique \
      --soloCellFilter EmptyDrops_CR \
      --soloUMIfiltering MultiGeneUMI_CR \
      --outSAMunmapped Within \
      --outFilterMismatchNmax 20 \
      --outFilterMismatchNoverLmax 0.1 \
      --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAA \
      --outFileNamePrefix $outdir/3prime_alignment/ \
      --outSAMmultNmax 1
'

#Align Internal reads as bulk

#internal_R1_barcoded=${internal_R1%.fastq.gz}.barcoded.fastq.gz
#internal_R2_barcoded=${internal_R2%.fastq.gz}.barcoded.fastq.gz

internal_R1_barcoded=$input_dir_internal/internal_barcoded_R1.fastq.gz
internal_R2_barcoded=$input_dir_internal/internal_barcoded_R2.fastq.gz
: '
STAR \
        --genomeDir $reference \
	--limitBAMsortRAM 50000000000 \
        --runThreadN 42 \
	--outSAMattrRGline ID:internal \
        --readFilesIn $internal_R1_barcoded $internal_R2_barcoded \
        --readFilesCommand zcat \
        --outSAMattributes GX GN \
        --outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within KeepPairs \
        --outFilterMismatchNmax 20 \
        --outFilterMismatchNoverLmax 0.1 \
	--clip3pAdapterSeq TGAAGCGGCGCACGAAAAACGCGAAAGCGTTTCAC TGGGCCGGTGCAGTTAATGTAGGGAAAGAGTGT \
	--clip3pAdapterMMp 0.1 0.1 \
	--outFileNamePrefix $outdir/internal_alignment/ \
	--outSAMmultNmax 1
'

#Align 5prime reads as bulk

#fiveprime_R1_barcoded=${fiveprime_R1%.fastq.gz}.barcoded.fastq.gz
#fiveprime_R2_barcoded=${fiveprime_R2%.fastq.gz}.barcoded.fastq.gz

fiveprime_R1_barcoded=$input_dir_5prime/fiveprime_barcoded_R1.fastq.gz
fiveprime_R2_barcoded=$input_dir_5prime/fiveprime_barcoded_R2.fastq.gz
: '
STAR \
        --genomeDir $reference \
	--genomeLoad LoadAndRemove \
	--limitBAMsortRAM 50000000000 \
        --runThreadN 42 \
	--outSAMattrRGline ID:5prime \
        --readFilesIn $fiveprime_R1_barcoded $fiveprime_R2_barcoded \
        --readFilesCommand zcat \
        --outSAMattributes GX GN \
        --outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within KeepPairs \
        --outFilterMismatchNmax 20 \
        --outFilterMismatchNoverLmax 0.1 \
        --clip3pAdapterSeq AAGCAGTGGTATCAACGCAGAGTACATGGG TGAAGCGGCGCACGAAAAACGCGAAAGCGTTTCAC \
        --clip3pAdapterMMp 0.1 0.1 \
        --outFileNamePrefix $outdir/5prime_alignment/ \
	--outSAMmultNmax 1
'


##merge BAM files and attach RG tag

mkdir -p $outdir/merged/
: '
samtools merge $outdir/merged/merged.bam \
	-f \
	-c \
	-p \
	-@ 32 \
	$outdir/3prime_alignment/Aligned.sortedByCoord.out.bam \
	$outdir/internal_alignment/Aligned.sortedByCoord.out.bam \
	$outdir/5prime_alignment/Aligned.sortedByCoord.out.bam
'



##add missing Gene alignents 
TMPDIR=$outdir/merged/_tmp

mkdir -p $TMPDIR

export TMPDIR

python /home/derek/analysis_10/SLR/scSLR/scripts/add_gene_tag_parallel.py \
	$GTF_FILE \
	$outdir/merged/merged.bam \
	12


##copy barcodes from read names to SAM tag 

samtools view -h $outdir/merged/merged.bam | awk '{if($1 ~ /^@/) print $0; else { for(i=1;i<=NF;i++){BAR=substr($1,1,18)} print $0,"TS:Z:"BAR}}' | samtools view -Sb - > $outdir/merged/merged.gene_tagged.sinto.bam



: '

samtools index \
	$outdir/merged/merged.gene_tagged.sinto.bam \
	-@ 32


## identify bead barcode groups using UMItools group

umi_tools group \
	--per-gene \
	--gene-tag GX \
	--stdin $outdir/merged/./marged.gene_tagged.sinto.3prime.firstInPair.sort.bam \
	--log $outdir/merged/group_log.txt \
	--method unique \
	--extract-umi-method tag \
	--umi-tag TS \
	--output-bam \
	--stdout $outdir/merged/marged.gene_tagged.sinto.3prime.firstInPair.sort.UMItooled.bam \
	--skip-tags-regex "[-]" \
	--no-sort-output \
	--timeit $outdir/merged/umi_tools_time.txt \
#	--unmapped-reads

samtools index \
	$outdir/merged/marged.gene_tagged.sinto.3prime.firstInPair.sort.UMItooled.bam \
	-@ 32


python share_UMIs.py \
	$outdir/merged/marged.gene_tagged.sinto.3prime.firstInPair.sort.UMItooled.sort.bam
'
echo "SLR pipeline completed"



