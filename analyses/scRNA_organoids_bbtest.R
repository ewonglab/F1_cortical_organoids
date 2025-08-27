

# Loading libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(Matrix)
  library(Matrix.utils)
  library(readxl)
  library(BiocParallel)
  library(openxlsx)
  library(DESeq2)
})  


# The beta-binomial test is performed on each strain separately

# CASTxBl6
strain <- "CAST"
dir <- "/data/"
allelic_count <- readRDS(paste0(dir, strain, ".organoids.allelic.counts.exons.rds"))
annot <- read.xlsx("/data/To_share_Good_cells_ALL.xlsx", sheet = 1, colNames = T)
colData(allelic_count)$replicate <- gsub("_.*", "", colnames(allelic_count))
colData(allelic_count)$cell_type <- annot$cell_type[match(gsub(".*_", "", rownames(colData(allelic_count))), annot$cell_barcode_clean)]
cell_info <- as.data.frame(colData(allelic_count))
cell_list <- split(cell_info, f = cell_info$cell_type)

allelic_count_split <- list()
for (i in 1:length(cell_list)){
  allelic_count_split[[i]] <- allelic_count[,rownames(cell_list[[i]])]
}

proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)
a1 <- list()
a2 <- list()
a1_norm <- list()
a2_norm <- list()
dat <- list()
for (i in 1:length(allelic_count_split)){
    a1[[i]] <- aggregate.Matrix(t(assays(allelic_count_split[[i]])[['a1']]),
                                groupings = colData(allelic_count_split[[i]])$replicate, fun = "sum")
    a2[[i]] <- aggregate.Matrix(t(assays(allelic_count_split[[i]])[['a2']]),
                                groupings = colData(allelic_count_split[[i]])$replicate, fun = "sum")
    a1[[i]] <- as.matrix(t(a1[[i]]))
    a2[[i]] <- as.matrix(t(a2[[i]]))
    colnames(a1[[i]]) <- paste0("BL6_", proper(strains[i]),"xBl6_", colnames(a1[[i]]))
    colnames(a2[[i]]) <- paste0(strains[i], "_", proper(strains[i]),"xBl6_", colnames(a2[[i]]))
    count_mat[[i]] <- aggregate.Matrix(t(assays(allelic_count_split[[i]])[['tot']]), 
                                       groupings = olData(allelic_count_split[[i]])$replicate, fun = "sum")
    count_mat[[i]] <- as.matrix(t(count_mat[[i]]))
    df[[i]] <- data.frame(clone = colnames(count_mat[[i]]))
    dds[[i]] <- DESeqDataSetFromMatrix(countData = count_mat[[i]], colData = df[[i]], design = ~ clone)
    dds[[i]] <- estimateSizeFactors(dds[[i]])
    sizeFactor[[i]] <- dds[[i]]$sizeFactor
    #normalising the counts
    a1_norm[[i]] <- sweep(a1[[i]], 2, sizeFactor, FUN = "/")
    a2_norm[[i]] <- sweep(a2[[i]], 2, sizeFactor, FUN = "/")
    colnames(a1_norm[[i]]) <- paste0("BL6_", proper(strain),"xBl6_", colnames(a1_norm[[i]]))
    colnames(a2_norm[[i]]) <- paste0(strain, "_", proper(strain),"xBl6_", colnames(a2_norm[[i]]))
    dat[[i]] <- cbind(a1_norm[[i]], a2_norm[[i]])

}

#Gliogenic progenitor cells do not appear after the treatment in clones 1 and 2, therefore those counts are removed
dat[[2]] <- dat[[2]][,-c(1:2,5:6)]

#generate total allelic count and allelic ratio tables
#sum first column with fifth, second with sixth, etc.
tot_counts_ct <- list()
genes_select <- list()
a1_counts_ct <- list()
ai_ratio_ct <- list()
for (i in 1:length(allelic_count_split)){
  tot_counts_ct[[i]] <- sapply(1:(ncol(dat[[i]])/2), function(x){rowSums(dat[[i]][,c(x,x+ncol(dat[[i]])/2)])})
  colnames(tot_counts_ct[[i]]) <- gsub("^[^_]*_", "", colnames(dat[[i]])[1:ncol(tot_counts_ct[[i]])])
  #select genes that have at least 10 counts across all replicates
  genes_select[[i]] <- rownames(tot_counts_ct[[i]])[rowSums(tot_counts_ct[[i]] >= 10) == ncol(tot_counts_ct[[i]])]
  tot_counts_ct[[i]] <- tot_counts_ct[[i]][genes_select[[i]],]
  a1_counts_ct[[i]] <- dat[[i]][genes_select[[i]], 1:ncol(tot_counts_ct[[i]])]
  ai_ratio_ct[[i]] <- a1_counts_ct[[i]]/tot_counts_ct[[i]]
  colnames(ai_ratio_ct[[i]]) <- gsub("^[^_]*_", "", colnames(dat[[i]])[1:ncol(tot_counts_ct[[i]])])
}


# Loading a list sex chromosome and imprinted genes
genesXY <- read.table(XY_genes)
genesIMPR <- read.xlsx(imprinted_genes)

#estimating global parameters of Beta-Binomial distribution
global_estim <- list()
for (i in 1:length(allelic_count_split)){
  global_estim[[i]] <- glob_disp(a1_counts_ct[[i]], tot_counts_ct[[i]], genesXY, genesIMPR)
}

param_estims <- list()
param_reestim_ct <- list()
bb_test_res <- list()
for (i in 1:length(allelic_count_split)){
  param_estims[[i]] <- estim_params_bulk(a1_counts_ct[[i]], tot_counts_ct[[i]])
  #change delta=2000 for SPRETxBl6 strain
  param_reestim_ct[[i]] <- correct_theta(param_estims[[i]], delta = 1000)
  bb_test_res[[i]] <- beta_binom_test_bulk(a1_counts_ct[[i]], tot_counts_ct[[i]], param_reestim_ct[[i]], global_estim[[i]])
  #no shrinkage
  bb_test_res[[i]]$fdr_orig <- p.adjust(bb_test_res[[i]]$pval_orig, method = "fdr")
  # with shrinkage
  bb_test_res[[i]]$fdr_shrunk <- p.adjust(bb_test_res[[i]]$pval_adj, method = "fdr")
  bb_test_res[[i]] <- bb_test_res[[i]][order(bb_test_res[[i]]$fdr_shrunk),]
  names(bb_test_res)[i] <- names(cell_list)[i]
} 

names(bb_test_res)[5] <- c('Intermediate neuronal prog')
names(bb_test_res)[7] <- c('Oligodendrocyte prog')

res_dir <- "/bbtest_res/"
filename <- paste0(res_dir, strain, "_bbtest_res.xlsx")
excel <- createWorkbook(filename)
for (i in 1:length(bb_test_res)){
  addWorksheet(excel, names(bb_test_res)[i])
  writeData(excel, sheet = i, bb_test_res[[i]], rowNames=T, colNames=T)
  saveWorkbook(excel, file = filename, overwrite = TRUE)
}



