#!/bin/bash

IN_BAM=$1

root=${IN_BAM%.bam}

bamCoverage \
	-b $IN_BAM \
	-o $root.bw \
	-p 4 \
	--ignoreDuplicates


