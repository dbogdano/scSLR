#!/bin/bash

IN_BAM=$1
OUT_BAM=$2

SUFFIX=${IN_BAM%"$.bam"}

picard UmiAwareMarkDuplicatesWithMateCigar --INPUT $IN_BAM --METRICS_FILE output_duplicate_metrics.txt --OUTPUT $OUT_BAM --UMI_METRICS_FILE UMI.txt --IGNORE_MISSING_MATES true --BARCODE_TAG CB --UMI_TAG_NAME UB --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT
	

