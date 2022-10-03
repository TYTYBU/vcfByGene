#!/bin/bash -l

BED_PATH='bed'
BED_PATH_LOCAL='mane_gene_bed'
VCF_OUT_PATH='vcf_out'
VCF_PATH="Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ final\ release"
dx mkdir -p $BED_PATH
dx mkdir -p $VCF_OUT_PATH

for i in {0..20} # process block 0-20
do
	for f in $BED_PATH_LOCAL/c1_b${i}_*.bed # edit here to process different chromosome and blocks
  do
  	f2=$(basename $f .bed)
  	dx upload --brief --wait --path $BED_PATH/$f2.bed $f

  	IFS='_' read -ra ADDR <<< "$f2"
  	BLK_STR="${ADDR[0]}_${ADDR[1]}"
  	GENE_STR=${ADDR[2]}

  	VCF_PREFIX="ukb23157_${BLK_STR}_v1"
  	BED_FILE=$BED_PATH/${f2}.bed

  	CMD_STRING="bcftools view --threads 4 -O z -R "$(basename $BED_FILE)" ${VCF_PREFIX}.vcf.gz > ${GENE_STR}.vcf.gz"
  	echo $CMD_STRING

  	dx run swiss-army-knife -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz" -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz.tbi" -iin="$BED_FILE" -icmd="$CMD_STRING" --destination vcf_out -y
  done
done

# for f in $BED_PATH_LOCAL/c1_b1_*.bed
# do
# 	f2=$(basename $f .bed)
# 	dx upload --brief --wait --path $BED_PATH/$f2.bed $f
#
# 	IFS='_' read -ra ADDR <<< "$f2"
# 	BLK_STR="${ADDR[0]}_${ADDR[1]}"
# 	GENE_STR=${ADDR[2]}
#
# 	VCF_PREFIX="ukb23157_${BLK_STR}_v1"
# 	BED_FILE=$BED_PATH/${f2}.bed
#
# 	CMD_STRING="bcftools view --threads 4 -O z -R "$(basename $BED_FILE)" ${VCF_PREFIX}.vcf.gz > ${GENE_STR}.vcf.gz"
# 	echo $CMD_STRING
#
# 	dx run swiss-army-knife -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz" -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz.tbi" -iin="$BED_FILE" -icmd="$CMD_STRING" --destination vcf_out -y
# done
