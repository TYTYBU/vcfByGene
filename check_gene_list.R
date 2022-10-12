library(openxlsx)
library(tidyverse)

mane = read.table("mane_track_UCSC_v1.0.txt", header = F, sep = "\t")

# pool genes from Chris, Rich and Vineel
temp = read.table("gene_list_Rich.txt", header = F)
temp2 = read.table("gene_list_Chris.txt", header = F)
temp3 = read.table("gene_list_Vineel.txt", header = F)
temp4 = read.table("gene_list_Vineel_2.txt", header = F)
temp = rbind(temp, temp2, temp3, temp4)
temp = unique(unlist(temp))

df = data.frame(gene_name = temp, gene_name2 = NA, bed_source = "MANE")
ind = which(df$gene_name %in% mane$V19 == F)
df$bed_source[ind] = NA

# some genes not found in MANE have alt gene names
ind = which(df$gene_name %in% c("AKAP2", "HIST1H1C", "C10orf113", "ACPP"))
df$gene_name2[ind] = c("PALM2AKAP2", "H1-2", "NEBL", "ACP3")
df$bed_source[ind] = "MANE"

temp = df$gene_name
temp[ind] = c("PALM2AKAP2", "H1-2", "NEBL", "ACP3")
write.table(temp, "genes_on_MANE.txt", quote = F, row.names = F, col.names = F)
write.xlsx(df, file = "selected_gene_summary.xlsx")


# gene list from Artur
temp = read_csv("resources/genebass_all-genes_Artur.csv")

df = data.frame(gene_name = temp$Gene, gene_name2 = NA, bed_source = "MANE")
ind = which(df$gene_name %in% mane$V19 == F)

