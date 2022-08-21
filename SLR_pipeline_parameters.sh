#!/bin/bash


### paramters for SLR pipeline to be source


output_directory=/media/chang/HDD-7/derek//SLR/GeneFull__6/GW23_1/

software=/home/derek/bin/



input_dir_3prime=/media/chang/HDD-7/derek/SLR/GeneFull__6/GW23_1/3prime

input_dir_internal=/media/chang/HDD-7/derek/SLR/GeneFull__6/GW23_1/internal

input_dir_5prime=/media/chang/HDD-7/derek/SLR/GeneFull__6/GW23_1/5prime

#STAR=/home/derek/reference_cellranger/cellranger-7.0.0/lib/bin/STAR



##optimized reference
#reference=/media/chang/HDD-1/reference_cellranger/STARSolo-GRCh38-v40/
reference=/home/derek/analysis_2/reference/human_hp3_reference/star


##standard 10X hg38
#reference=/media/chang/HDD-1/reference_cellranger/STAR_GRCh38-2020-A/

#GTF_FILE=/media/chang/HDD-1/reference_cellranger/refdata-gex-GRCh38-v40/genes/genes.gtf
GTF_FILE=/home/derek/analysis_2/reference/human_hp3_reference/genes/genes.gtf




##For visium use visium whitelist:
#BC_whitelist=/media/chang/HDD-1/reference_cellranger/cellranger-6.1.2/lib/python/cellranger/barcodes/visium-v1.txt



BC_whitelist=/media/chang/HDD-1/reference_cellranger/cellranger-6.1.2/lib/python/cellranger/barcodes/3M-february-2018.txt

TellREAD_whitelist=/home/derek/analysis_7/SLR/TellSeq_whitelist






