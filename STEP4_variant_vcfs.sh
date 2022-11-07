#!/bin/bash -l

VCF_PATH="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_vcf_out"
BED_PATH="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_diff_regions"
REF_FILE="project-GGy3Bb0JqBj7zfxY8v4by61X:/GRCh38_resources/GRCh38_reference_genome.fa"
REF_FILE_INDEX="project-GGy3Bb0JqBj7zfxY8v4by61X:/GRCh38_resources/GRCh38_reference_genome.fa.fai"

NO_DIFF_GENES="CAV3 CSRP3 HEPH H1-2 MYL3 TTR"
GENE_VCFS=$(dx ls selected_genes_vcf_out)

for gene_file in $GENE_VCFS; do
  gene=${gene_file%.vcf.gz}

  gene_path="$VCF_PATH/$gene_file"
  bed_path="$BED_PATH/$gene.bed"

  cmd_str="bcftools view -T ^$gene.bed -Ou $gene.vcf.gz | bcftools norm -m - -f GRCh38_reference_genome.fa -Ou | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou | bcftools view --threads 4 -i 'MAF<=0.001 && MAC >=1 && F_MISSING<0.1 && F_PASS(DP>=10 & GT!=\"mis\")> 0.9' -Ou | bcftools view -G -Ou | bcftools annotate --threads 4 -x INFO -Oz > ${gene}_variants.vcf.gz"
  #echo $cmd_str
  dx run swiss-army-knife -iin="$gene_path" -iin="$REF_FILE" -iin="$REF_FILE_INDEX" -iin="$bed_path" -icmd="$cmd_str" --destination "selected_genes_filtered_variants" --tag "var_pt_bcftools" -y

done

for gene in $NO_DIFF_GENES; do

  gene_path="$VCF_PATH/$gene.vcf.gz"

  cmd_str="bcftools norm -m - -f GRCh38_reference_genome.fa -Ou $gene.vcf.gz | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou | bcftools view --threads 4 -i 'MAF<=0.001 && MAC >=1 && F_MISSING<0.1 && F_PASS(DP>=10 & GT!=\"mis\")> 0.9' -Ou | bcftools view -G -Ou | bcftools annotate --threads 4 -x INFO -Oz > ${gene}_variants.vcf.gz"
  #echo $cmd_str
  dx run swiss-army-knife -iin="$gene_path" -iin="$REF_FILE" -iin="$REF_FILE_INDEX" -icmd="$cmd_str" --destination "selected_genes_filtered_variants" --tag "var_pt_bcftools" -y

done

CMD_STR="bcftools view -H LDLR_variants.vcf.gz > LDLR.txt"
dx run swiss-army-knife -iin="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_filtered_variants/LDLR_variants.vcf.gz" -icmd="$CMD_STR" --destination "project-GGy3Bb0JqBj7zfxY8v4by61X:/hail_pipelines/LDLR_test"

# CMD="bcftools view -H A1CF_variants.vcf.gz | head -n 20 > A1CF_variants.txt"
# dx run swiss-army-knife -iin="project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_filtered_variants/A1CF_variants.vcf.gz" -icmd="$CMD" --destination "project-GGy3Bb0JqBj7zfxY8v4by61X:/selected_genes_filtered_variants"
