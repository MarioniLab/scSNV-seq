# This script counts barcode occurences from sequences, using read groups and 
# sequences extracted from reads mapped to the gRNA and iBAR amplicons by
# the MissionBio tapestri pipeline

sample_name=commandArgs(trailingOnly=TRUE)
sequence_length<- 151
library(Biostrings)
library(dplyr)

data_all <- read.table(paste0("tapestri_",sample_name,".sequences.txt"))
sequences <- data_all$V1
cells_ind_1 <- grepl("RG", data_all$V2)
cells_ind_2 <- grepl("RG", data_all$V3)
cells <- data_all$V2
cells[cells_ind_2] <- data_all$V3[grepl("RG", data_all$V3)]
cells <- substring(cells,6,25)
df_cells_sequences <- data.frame(cell=cells,sequence=sequences)
yy <- width(sequences) == sequence_length
df_cells_sequences <- df_cells_sequences[yy,]

df_cells_sequences_puro <- df_cells_sequences[agrepl(toupper("gcagctgaagcttaccatga"), df_cells_sequences$sequence, max.distance = 1),]

df_cells_sequences_gRNA <- df_cells_sequences[agrepl(toupper("gtttaagagctatgctgga"), df_cells_sequences$sequence, max.distance = 1),]

find_puro <-  vmatchPattern(toupper("gcagctgaagcttaccatga"), df_cells_sequences_puro$sequence, max.mismatch=1)
puro_ends <- unlist(find_puro@ends)
df_cells_sequences_puro$puro_barcode <- subseq(df_cells_sequences_puro$sequence,puro_ends-31,puro_ends-20)

df_puro_count <- df_cells_sequences_puro[c("puro_barcode","cell")] %>% group_by_all() %>% count()
saveRDS(df_puro_count,file= paste0("MB_puroBC_",sample_name,".rds"))

find_gRNA <-  vmatchPattern(toupper("gtttaagagctatgctgga"), df_cells_sequences_gRNA$sequence, max.mismatch=1)
gRNA_ends <- unlist(find_gRNA@ends)
xx <- gRNA_ends > 90
xx[gRNA_ends > 95] <- FALSE
df_cells_sequences_gRNA <- df_cells_sequences_gRNA[xx,]
df_cells_sequences_gRNA$gRNA <- subseq(df_cells_sequences_gRNA$sequence,gRNA_ends[xx]-37,gRNA_ends[xx]-19)
df_cells_sequences_gRNA$iBAR <- subseq(df_cells_sequences_gRNA$sequence,gRNA_ends[xx]+1,gRNA_ends[xx]+6)

df_gRNA_iBAR_count <- df_cells_sequences_gRNA[c("gRNA","iBAR","cell")] %>% group_by_all() %>% count()
saveRDS(df_gRNA_iBAR_count,file= paste0("MB_gRNA_iBAR_",sample_name,".rds"))

df_gRNA_iBAR_count_no_cells <- df_cells_sequences_gRNA[c("gRNA","iBAR")] %>% group_by_all() %>% count()
saveRDS(df_gRNA_iBAR_count_no_cells,file= paste0("MB_gRNA_iBAR_mapping_",sample_name,".rds"))

df_gRNA_count <- df_cells_sequences_gRNA[c("gRNA","cell")] %>% group_by_all() %>% count()
saveRDS(df_gRNA_count,file= paste0("MB_gRNA_",sample_name,".rds"))

