#!/bin/bash -l

VCF_PATH="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_vcf_out"
BED_PATH="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_diff_regions"
REF_FILE="project-GGy3Bb0JqBj7zfxY8v4by61X:/GRCh38_resources/GRCh38_reference_genome.fa"
REF_FILE_INDEX="project-GGy3Bb0JqBj7zfxY8v4by61X:/GRCh38_resources/GRCh38_reference_genome.fa.fai"

GENE_VCFS=$(dx ls selected_genes_vcf_out)

for gene_file in $GENE_VCFS; do
  gene=${gene_file%.vcf.gz}

  gene_path="$VCF_PATH/$gene_file"
  bed_path="$BED_PATH/$gene.bed" #note: fails for genes where this file is empty!

  cmd_str="bcftools view -T ^$gene.bed -Ou $gene.vcf.gz | bcftools norm -m - -f GRCh38_reference_genome.fa -Ou | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou | bcftools view --threads 4 -i 'MAF<=0.001 && MAC >=1 && F_MISSING<0.1 && F_PASS(DP>=10 & GT!=\"mis\")> 0.9' -Ou | bcftools query -i 'GT=\"RA\"|GT=\"AR\"|GT=\"AA\"' -f '%CHROM  %POS %REF %ALT [%SAMPLE|]\n' > $gene.ssv"
  #echo $cmd_str
  dx run swiss-army-knife -iin="$gene_path" -iin="$REF_FILE" -iin="$REF_FILE_INDEX" -iin="$bed_path" -icmd="$cmd_str" --destination "selected_genes_variant_patient" --tag "var_pt_bcftools" -y

done

# MAP_PATH="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_variant_patient"
# VARIANT_MAPS=$(dx ls selected_genes_variant_patient)
#
# for patient_file in $VARIANT_MAPS; do
#   dx run gzip -ifile="$MAP_PATH/$patient_file" -icompression_level=9
# done
