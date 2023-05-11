
library(Matrix)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(mixsmsn)
library(biomaRt)

read_in_filtered_matrix <- function(name_path)
{# for cellranger output
  matrix_dir = paste(name_path,"/filtered_feature_bc_matrix/",sep="")
  barcode_path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features_path <- paste0(matrix_dir, "features.tsv.gz")
  matrix_path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat_filtered <- readMM(file = matrix_path)
  feature_names = read.delim(features_path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode_names = read.delim(barcode_path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat_filtered) = barcode_names$V1
  rownames(mat_filtered) = feature_names$V1
  return(mat_filtered)
}

.concatenate_vector_of_strings <- function(v){
  if (length(v) <= 1){
    return(v)
  }else{
    w <- v[1]
    for (j in 2:length(v)){
      w <- paste0(w,"-",v[j])
    }
    return(w)
  }
}

.assign_barcodes <- function(barcode_count_matrix,thresh_lower,thresh_upper){
  #keep only cells without barcodes with (umi) counts between the lower and upper bounds
  ind1 <- apply(barcode_count_matrix,2,function(x) return(sum(x>thresh_upper) == sum(x>thresh_lower)))
  temp <- barcode_count_matrix[,ind1]
  barcode_list <- apply(temp,2,function(x) rownames(temp)[x > thresh_upper])
  output <- c()
  output$barcode_assignment <- lapply(barcode_list,.concatenate_vector_of_strings)
  output$barcode_assignment <- unlist(output$barcode_assignment[unlist(lapply(output$barcode_assignment, function(x) !(is.null(x))))])
  output$barcode_assignment_unique <- NA
  output$barcode_assignment_unique[!(grepl("-", output$barcode_assignment))] <- output$barcode_assignment[!(grepl("-", output$barcode_assignment))] 
  output$barcode_distribution <- table(unlist(barcode_list))
  output$barcode_distribution_unique <- table(output$barcode_assignment_unique)
  return(output)
}


#' @title barcode_calling
#' @param barcode_matrix b x c matrix of counts, where c is number of cells and b the number of barcodes, the colnames and rownames of the matrix are the cells and barcodes respectively,
#' the element in row j and column k of the matrix is the number of (UMI) counts for barcode j in cell k 
#' @param sample_name name of the sample
#' @param barcode_name name of the barcode
#' @param thresh minimum counts of UMI count threshold for potential barcode assignment,
#' default value 1; the value it should be set to depends on the sequencing depth, 
#' set to a value that can safely be assumed to be ambient noise, but also high enough
#' in case sequencing is deeper to avoid inefficiency and modelling errors by huge numbers
#' of low ambient counts included in the mixture modelling
#' @param g number of components of the mixture model. Generally 3 components are recommended, as the 
#' 'signal component' corresponding to gRNAs assigned to cells tends to be bimodal. For specific data sets
#' g=2 or other values may improve assignment, but changing g should be based on prior inspection of the data. 
#' @export
barcode_calling <- function(barcode_matrix,sample_name="",barcode_name="gRNA",thresh=1,g=3){
  temp <- barcode_matrix[barcode_matrix  > thresh]
  temp <- temp[!(is.na(temp))]
  temp <- sort(temp)
  fit <- smsn.mix(log2(temp), nu = 0, g = g, get.init = TRUE, criteria = TRUE,group = TRUE, family = "Skew.normal", calc.im=FALSE,obs.prob=TRUE)
  
  print(mix.hist(log2(temp),fit))
  ind1 <- which.min(fit$mu)
  thresh_lower <- temp[max(which(fit$obs.prob[,ind1] > 0.9 ))]
  if (thresh_lower < thresh){
    thresh_lower <- thresh
  }
  thresh_upper <- temp[min(which(fit$obs.prob[,ind1] <  0.1))]
  df_assign <- data.frame(prob_cluster1 = fit$obs.prob[,ind1],log_barcode_counts = log2(temp))
  df_thresh = data.frame(xintercept = c(log2(thresh_lower), log2(thresh_upper)), thresholds = c("thresh_lower", "thresh_upper"))
  
  p1 <- ggplot(df_assign,aes(x=log_barcode_counts,y=prob_cluster1) )+ geom_point(size=3,alpha=0.01) +
    geom_vline(aes(xintercept = xintercept, color = thresholds), df_thresh,size=3,alpha=0.5) + xlab("") + ylab("probability of being\nin left cluster")+
    theme_classic(base_size=13)+scale_color_manual(values=c("thresh_lower"="darkblue","thresh_upper"="darkred"))+ 
    theme(legend.key = element_rect(fill = "white", colour = "black"))+ theme(legend.position="bottom")+
    ggtitle(paste0("barcode ",sample_name))
  
  
  barcode_assignment <-.assign_barcodes(barcode_matrix,thresh_lower,thresh_upper)
  saveRDS(barcode_assignment,file=paste0("barcode_assignment_",sample_name,".rds"))
  
  
  df_assign <- data.frame(log_barcode_counts = log2(temp))
  df_thresh = data.frame(xintercept = c(log2(thresh_lower), log2(thresh_upper)), thresholds = c("thresh_lower", "thresh_upper"))
  
  p3 <- ggplot(df_assign,aes(x=log_barcode_counts))+ geom_histogram() +
    geom_vline(aes(xintercept = xintercept, color = thresholds), df_thresh,size=2,alpha=0.5) + ylab("")+
    theme_classic(base_size=13)+scale_color_manual(values=c("thresh_upper"="darkred","thresh_lower" = "darkblue"))+ 
    theme(legend.key = element_rect(fill = "white", colour = "black"))+ theme(legend.position="bottom")
  nr_barcodes <- sapply(barcode_assignment$barcode_assignment,function(x) sum(grepl("-",x))+1)
  p4 <- ggplot(mapping=aes(x=as.factor(nr_barcodes))) + geom_histogram(stat="count") + theme_classic(base_size=13)  +
    xlab(paste0("number of ",barcode_name," per cell")) + ylab("number of cells")
  
  
  pp <- (p1+p3+p4)+ plot_annotation(title=sample_name)
  print(pp)
  output <- c()
  output$thresh_lower <- thresh_lower
  output$thresh_upper <- thresh_upper
  output$barcode_assignment <- barcode_assignment
  return(output)
}

cell_cycle_scoring <- function(sce)
{
  sce$S.Score <- rep(NA,ncol(sce))
  sce$G2M.Score <- rep(NA,ncol(sce))
  sce$phase <- rep(NA,ncol(sce))
  s.genes = cc.genes$s.genes
  g2m.genes = cc.genes$g2m.genes
  SeuratObj <- Seurat::CreateSeuratObject(as.matrix(counts(sce)))
  SeuratObj <- Seurat::NormalizeData(SeuratObj)
  SeuratObj <- Seurat::FindVariableFeatures(SeuratObj, selection.method = "vst")
  SeuratObj <- Seurat::ScaleData(SeuratObj, features = rownames(SeuratObj))
  SeuratObj <- Seurat::CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  sce$S.Score <- SeuratObj$S.Score
  sce$G2M.Score <- SeuratObj$G2M.Score
  sce$phase  <- SeuratObj$Phase
  return(sce)
}


#' @title QC_mRNA_outlier
#' QC for mRNA modality, removes low outliers for the total count, low outliers for the number of detected features 
#' and high outliers for the percentage of counts from mitochondrial genes 
#' @param mRNA_counts matrix of UMI counts (rows:genes, columns:cells), rownames need to be the ensembl gene ids
#' @export
QC_mRNA_outlier <- function(mRNA_counts,file_name="")
{
  sce <- SingleCellExperiment(
    assays = list(counts = mRNA_counts), colData = colnames(mRNA_counts))
  temp <- query(AnnotationHub(), "EnsDb.Hsapiens.")
  edb <- AnnotationHub()[["AH95744"]]
  chr_loc <- mapIds(edb, keys=rownames(sce),
                    keytype="GENEID", column="SEQNAME")
  is_mito <- which(chr_loc=="MT")
  df <- perCellQCMetrics(sce, subsets=list(Mito=is_mito))
  A <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
  colSums(as.matrix(A))
  discard <- A$discard
  colData(sce) <- cbind(colData(sce), df)
  sce$discard <- discard
  plot_list <- list()
  plot_list[[1]] <- plotColData(sce,  y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count")
  plot_list[[2]] <- plotColData(sce, y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features")
  plot_list[[3]] <- plotColData(sce,  y="subsets_Mito_percent",
                                colour_by="discard") + ggtitle("Mito percent")
  p <- marrangeGrob(plot_list,nrow=3,ncol=1)
  print(p)
  ggsave(p,filename=paste0(file_name,"QC_plot.pdf"))
  return(discard)
}


#' @title getHVGs
#' Identifies highly variable genes in the data.
#' @param sce SingleCellExperiments with logcounts assay
#' @param min_mean minimum mean expression cutoff for a gene to be 
#' potentially included in the group of highly variable genes
#' @param FDR FDR cutoff for calling highly variable genes
#' @export
getHVGs <- function(sce,min_mean = 0.001,FDR = 0.01)
{
  stats <- modelGeneVar(sce)
  stats <- stats[stats$mean > min_mean,]
  top_genes <- rownames(stats[stats$FDR < FDR,])
  return(top_genes)
}

read_in_raw_matrix <- function(name_path)
{# for cellranger output
  matrix_dir = paste(name_path,"/raw_feature_bc_matrix/",sep="")
  barcode_path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features_path <- paste0(matrix_dir, "features.tsv.gz")
  matrix_path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat_raw <- readMM(file = matrix_path)
  feature_names = read.delim(features_path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode_names = read.delim(barcode_path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat_raw) = barcode_names$V1
  rownames(mat_raw) = feature_names$V1
  return(mat_raw)
}

#' @title getmode
#' Helper function to transfer UMAP coordinates and cluster labels 
#' to a mapped data set. 
getmode <- function(v, dist) {
  # taken from
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

