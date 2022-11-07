#!/bin/bash -l
BED_PATH='Artur_genes_bed'
BED_PATH_LOCAL='temp' # need to be the same as the mane_bed.tar.gz extraction folder
VCF_OUT_PATH='Artur_genes_vcf_out'
VCF_PATH="Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ final\ release"
dx mkdir -p $BED_PATH
dx mkdir -p $VCF_OUT_PATH

tar -zxf mane_bed.tar.gz # extract bed files to ./temp/

# process VCFs by a list of genes
RESOURCE_PATH='resources'

# GENE_LIST=$RESOURCE_PATH/genes_not_on_MANE.txt
# while read gene; do
# 	for f in $BED_PATH_LOCAL/*_${gene}_*.bed # edit here to process different chromosome and blocks
#   do
# 		if [ -f "$f" ]; then
# 			# echo $f
# 	  	f2=$(basename $f .bed)
# 	  	dx upload --brief --wait --path $BED_PATH/$f2.bed $f
#
# 	  	IFS='_' read -ra ADDR <<< "$f2"
# 	  	BLK_STR="${ADDR[0]}_${ADDR[1]}"
# 	  	GENE_STR=${ADDR[2]}_${ADDR[3]}_${ADDR[4]} # refseq_id after gene symbol
#
# 	  	VCF_PREFIX="ukb23157_${BLK_STR}_v1"
# 	  	BED_FILE=$BED_PATH/${f2}.bed
#
# 	  	CMD_STRING="bcftools view --threads 4 -O z -R "$(basename $BED_FILE)" ${VCF_PREFIX}.vcf.gz > ${GENE_STR}.vcf.gz"
# 	  	echo $CMD_STRING
#
# 	  	dx run swiss-army-knife -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz" -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz.tbi" -iin="$BED_FILE" -icmd="$CMD_STRING" --destination $VCF_OUT_PATH -y
# 		else
# 			echo "$f does not exist."
# 		fi
#   done
# done <${GENE_LIST}

GENE_LIST=$RESOURCE_PATH/genes_on_MANE_Artur.txt
while read gene; do
	for f in $BED_PATH_LOCAL/*_${gene}.bed # edit here to process different chromosome and blocks
  do
		if [ -f "$f" ]; then
			# echo $f
	  	f2=$(basename $f .bed)
	  	dx upload --brief --wait --path $BED_PATH/$f2.bed $f

	  	IFS='_' read -ra ADDR <<< "$f2"
	  	BLK_STR="${ADDR[0]}_${ADDR[1]}"
	  	GENE_STR=${ADDR[2]}

	  	VCF_PREFIX="ukb23157_${BLK_STR}_v1"
	  	BED_FILE=$BED_PATH/${f2}.bed

	  	CMD_STRING="bcftools view --threads 4 -O z -R "$(basename $BED_FILE)" ${VCF_PREFIX}.vcf.gz > ${GENE_STR}.vcf.gz"
	  	echo $CMD_STRING

	  	dx run swiss-army-knife -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz" -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz.tbi" -iin="$BED_FILE" -icmd="$CMD_STRING" --destination $VCF_OUT_PATH -y
		else
			echo "$f does not exist."
		fi
  done
done <${GENE_LIST}



# # process VCFs by blocks
# for i in {0..20} # process block 0-20
# do
# 	for f in $BED_PATH_LOCAL/c1_b${i}_*.bed # edit here to process different chromosome and blocks
#   do
#   	f2=$(basename $f .bed)
#   	dx upload --brief --wait --path $BED_PATH/$f2.bed $f
#
#   	IFS='_' read -ra ADDR <<< "$f2"
#   	BLK_STR="${ADDR[0]}_${ADDR[1]}"
#   	GENE_STR=${ADDR[2]}
#
#   	VCF_PREFIX="ukb23157_${BLK_STR}_v1"
#   	BED_FILE=$BED_PATH/${f2}.bed
#
#   	CMD_STRING="bcftools view --threads 4 -O z -R "$(basename $BED_FILE)" ${VCF_PREFIX}.vcf.gz > ${GENE_STR}.vcf.gz"
#   	echo $CMD_STRING
#
#   	dx run swiss-army-knife -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz" -iin="$VCF_PATH/${VCF_PREFIX}.vcf.gz.tbi" -iin="$BED_FILE" -icmd="$CMD_STRING" --destination vcf_out -y
#   done
# done

rm -rf $BED_PATH_LOCAL
rm -rf mane_bed.tar.gz
