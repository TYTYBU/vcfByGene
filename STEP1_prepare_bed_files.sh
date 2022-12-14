#!/bin/bash -l

BED_PATH_LOCAL='temp'
RESOURCE_PATH='resources'
mkdir -p $BED_PATH_LOCAL
perl mane2bed.pl $RESOURCE_PATH/mane_track_UCSC_v1.0.txt $RESOURCE_PATH/pvcf_blocks.txt $BED_PATH_LOCAL
perl refseq2bed.pl $RESOURCE_PATH/refseq_7genes.txt $RESOURCE_PATH/pvcf_blocks.txt $BED_PATH_LOCAL
# tar czf mane_bed.tar.gz $BED_PATH_LOCAL/*.bed
# rm -rf $BED_PATH_LOCAL
