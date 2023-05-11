#This script summarises per cell GATK vcf file output to a list in R


library(dplyr)

vcf_files <- list.files(getwd(),pattern="vcf")
vcf_files <- setdiff(vcf_files,vcf_files[grepl("tbi",vcf_files)])
genotype_per_cell <- list()

for (j in 1:length(vcf_files)){
  temp <- read.vcfR(cf_files[j]),verbose=FALSE)
  fix <- temp@fix
  gt <- temp@gt
  df_temp <- data.frame(chrom=fix[,1], pos=fix[,2],
                        ref=fix[,4],alt=fix[,5])
  df_temp$type <- sapply(gt[,2],function(x) return(substr(x,1,5)))
  ind_chr1 <- df_temp$chrom=="chr1"
  df_temp <- df_temp[ind_chr1,]
  df_fix <- as.data.frame(fix)
  df_fix <- df_fix[which(ind_chr1),]
  genotype_per_cell[[j]] <- df_temp
}

names_cells <- sapply(vcf_files,function(x) return(substring(x,1,20)))
names(genotype_per_cell) <- names_cells
saveRDS(genotype_per_cell,file="genotype_per_cell_all.rds")

