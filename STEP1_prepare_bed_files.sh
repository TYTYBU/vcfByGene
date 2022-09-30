#!/bin/bash -l

BED_PATH_LOCAL='mane_gene_bed'
RESOURCE_PATH='resources'
mkdir -p $BED_PATH_LOCAL

perl mane2bed.pl $RESOURCE_PATH/mane_track_UCSC_v1.0.txt $RESOURCE_PATH/pvcf_blocks.txt $BED_PATH_LOCAL
