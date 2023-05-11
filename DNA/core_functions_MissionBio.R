source("../core_functions.R")
library(dplyr)

#' @title subset_C_T
#' @param df_gt data frame with at least the following columns: 
#' ref-reference allele; alt-alternative allele
subset_C_T <- function(df_gt){
  xx <- df_gt$ref == "C" & df_gt$alt=="T"
  xx[df_gt$ref == "C" & grepl("T",df_gt$alt) & grepl(",",df_gt$alt) & nchar(df_gt$alt) <= 3] <- TRUE
  xx[df_gt$ref == "G" & df_gt$alt == "A"] <- TRUE
  xx[df_gt$ref == "G" &  grepl("A",df_gt$alt) & grepl(",",df_gt$alt) & nchar(df_gt$alt) <= 3] <- TRUE
  return(df_gt[as.logical(xx),])
}

  
#' @title assign_mutations_to_puro_barcodes
#' @param df_GT_min_cell_number data frame including cells for barcodes with a 
#' minimum number of cells per barcode; columns: GT-genotype, cell, and a column
#' containing the assigned barcode, specified in column_name_barcode
#' @param column_name_barcode name of the column of the data frame that specifies
#' the assigned barcode
#' @param cutoff_WT_nr maximum number of cells with each edited position for the barcode to 
#' be assigned WT
#' @param cutoff_WT_frequency maximum frequency for each edited position for the barcode to be assigned WT
#' @param cutoff_freq_pos minimum frequency of cells edited at a position for this position to be included
#' as part of the assigned genotype for the barcode
#' @param cutoff_freq_type minimum frequency of cells edited for a specific number of alleles to be included
#' as part of the assigned genotype for the barcode and specific position
#' @export
#' @examples 

assign_mutations_to_puro_barcodes <- function(df_GT_min_cell_number, column_name_barcode, cutoff_WT_nr=1, cutoff_WT_freq=0.1, 
                                              cutoff_nr_pos=2,cutoff_freq_pos=0.5,cutoff_freq_type=0.51){
  unique_barcode <- unique(df_GT_min_cell_number[[column_name_barcode]])
  GT_puro_BC <- rep("",length(unique_barcode))
  for (j in 1:length(unique_barcode)){
    temp_mutations <- unlist(sapply(df_GT_min_cell_number$GT[df_GT_min_cell_number[[column_name_barcode]] == unique_barcode[j]],function(x) strsplit(x,";")[[1]]))
    temp_positions <- sapply(temp_mutations,function(x) strsplit(x,"-")[[1]][1])
    temp_mutations_edited <- temp_mutations[temp_positions!=""]
    temp_positions <- temp_positions[temp_positions!=""]
    temp_AF <- sapply(temp_mutations_edited,function(x) strsplit(x,"-")[[1]][2])
    temp_positions_at_least_2_alleles <- temp_positions[temp_AF%in%c("0/1/1","1/1/1")]
    freq_positions <- table(temp_positions)/sum(df_GT_min_cell_number[[column_name_barcode]] == unique_barcode[j])
    if ((all(table(temp_positions)<= cutoff_WT_nr) | all(freq_positions <= cutoff_WT_freq))) {
      GT_puro_BC[j] <- "WT"
    }else{
      positions_keep <- names(table(temp_positions))[table(temp_positions) >= cutoff_nr_pos & freq_positions>=cutoff_freq_pos]
      positions_keep <- unique(positions_keep)
      mutations_keep <- c()
      if (length(positions_keep) > 0){
        for (k in 1:length(positions_keep)){
          mutations_temp1 <- temp_mutations[grepl(positions_keep[k],temp_mutations)]
          table_mutations_temp1 <- table(mutations_temp1)
          #freq_mutations_temp1 <- table_mutations_temp1/sum(df_GT_min_cell_number[[column_name_barcode]] == unique_barcode[j])
          freq_mutations_temp1 <- table_mutations_temp1/sum(table_mutations_temp1)
          if (max(table_mutations_temp1) >= cutoff_freq_type){
            mutations_keep <- c(mutations_keep,names(table_mutations_temp1)[table_mutations_temp1==max(table_mutations_temp1)])}
        }
        mutations_keep <- sort(mutations_keep)
        if (length(mutations_keep) >= 1){
          GT_puro_BC[j] <- mutations_keep[1]
          if (length(mutations_keep) >1){
            for (k in 2:length(mutations_keep)){
              GT_puro_BC[j] <- paste0(GT_puro_BC[j],";",mutations_keep[k])
            }
          }
        }
      }
    }
    
  }
  names(GT_puro_BC) <- unique_barcode
  return(GT_puro_BC)
}
