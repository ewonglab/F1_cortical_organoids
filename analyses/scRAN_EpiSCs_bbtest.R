

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


# Functions

#splitting STARSolo output by allele-specific counts
allele_split <- function(obj){
  
  a1 <- obj@assays[["RNA"]]@layers[["counts"]][!grepl("pseudo", rownames(obj@assays[["RNA"]])),]
  a2 <- obj@assays[["RNA"]]@layers[["counts"]][grepl("pseudo", rownames(obj@assays[["RNA"]])),]
  rownames(a1) <- rownames(obj@assays[["RNA"]])[!grepl("pseudo", rownames(obj@assays[["RNA"]]))]
  rownames(a2) <- rownames(obj@assays[["RNA"]])[grepl("pseudo", rownames(obj@assays[["RNA"]]))]
  colnames(a1) <- colnames(obj@assays[["RNA"]])
  colnames(a2) <- colnames(obj@assays[["RNA"]])
  rownames(a2) <- gsub("pseudo", "", rownames(a2))
  
  idx <- rownames(a2) %in% rownames(a1)
  genes <- rownames(a2)[idx]
  
  a1 <- a1[genes,]
  a2 <- a2[genes,]
  
  tot <- a1 + a2
  ai <- a1/tot
  
  SingleCellExperiment(assays = list(a1 = as.matrix(a1),
                                     a2 = as.matrix(a2),
                                     tot = as.matrix(tot),
                                     ai = as.matrix(ai)))
  
}

# Removing low expressed genes
filter_sc <- function(sce){
  filter <- rowSums(assays(sce)[['tot']] > 1) >= 10 
  sce <- sce[filter,]
}

# Converting SingleCellExperiment object to Seurat object
convert2Seurat <- function(obj){
  CreateSeuratObject(obj@assays@data@listData[["tot"]])
}

# preprocess data in Seurat
preprocess <- function(obj){
  obj  <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(obj, dims = 1:20, k.param = 15)
  obj <- FindClusters(obj, resolution = 0.5)
  obj <- RunUMAP(obj, dims = 1:20)
}


# Loading data
strains <- c("CAST", "MOLF", "PWK", "SPRET")
res_dir <- "/align/GRCm38/scRNA_EpiSC/"
filt_dir <- paste0(res_dir, strains, "/STARsoloIndex/", strains, "_STARsoloIndex_Solo.out/Gene/filtered")

# Loading STARSolo output

import_seurat <- function(dir){
  filt.data <- Read10X(dir)
  seur.data <- CreateSeuratObject(counts = filt.data, project = "EpiSC", min.cells = 5, min.features = 100)
}  

seur.obj <- sapply(filt_dir, import_seurat)
names(seur.obj) <- paste0(strains, "xBL6")

VlnPlot(seur.obj[[1]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + labs(caption = "CASTxBL6")
VlnPlot(seur.obj[[2]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + labs(caption = "MOLFxBL6")
VlnPlot(seur.obj[[3]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + labs(caption = "PWKxBL6")
VlnPlot(seur.obj[[4]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + labs(caption = "SPRETxBL6")


#Generating allelic counts object
allelic_counts <- sapply(seur.obj, allele_split)

#Remove lowly expressed genes which have less than 10 cells with any counts
allelic_count_filt <- sapply(allelic_counts, filter_sc)

# In the absence of real biological replicates, Seurat clusters will be used to pseudo-bulk the counts
seur.obj2 <- sapply(allelic_count_filt, convert2Seurat)

# To generate clusters data is pre-process using Seurat workflow
seur.obj2 <- sapply(seur.obj2, preprocess)

# Visualising clusters
DimPlot(seur.obj2[[1]], reduction = "umap") + labs(caption = "CASTxBl6")
DimPlot(seur.obj2[[2]], reduction = "umap") + labs(caption = "MOLFxBl6")
DimPlot(seur.obj2[[3]], reduction = "umap") + labs(caption = "PWKxBl6")
DimPlot(seur.obj2[[4]], reduction = "umap") + labs(caption = "SPRETxBl6")


# Adding cluster annotation to the SingleCellExperiment metadata
colData(allelic_count_filt[[1]])$seurat_clusters <- seur.obj2[[1]]@meta.data[["seurat_clusters"]]
colData(allelic_count_filt[[2]])$seurat_clusters <- seur.obj2[[2]]@meta.data[["seurat_clusters"]]
colData(allelic_count_filt[[3]])$seurat_clusters <- seur.obj2[[3]]@meta.data[["seurat_clusters"]]
colData(allelic_count_filt[[4]])$seurat_clusters <- seur.obj2[[4]]@meta.data[["seurat_clusters"]]


#Removing an outlier cluster with PWKxBl6 hybrid (90 cells in total)
allelic_count_filt[[3]] <- allelic_count_filt[[3]][,colData(allelic_count_filt[[3]])$seurat_clusters != "4"]


# Aggregate single-cell counts
groups <- list()
for (i in 1:length(allelic_count_filt)){
  groups[[i]] <- colData(allelic_count_filt[[i]])[, c("seurat_clusters")]
}

proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)
a1 <- list()
a2 <- list()
dat <- list()
for (i in 1:length(allelic_counts)){
  a1[[i]] <- aggregate.Matrix(t(assays(allelic_count_filt[[i]])[['a1']]), 
                              groupings = colData(allelic_count_filt[[i]])$seurat_clusters, fun = "sum")
  a2[[i]] <- aggregate.Matrix(t(assays(allelic_count_filt[[i]])[['a2']]), 
                              groupings = colData(allelic_count_filt[[i]])$seurat_clusters, fun = "sum")
  a1[[i]] <- as.matrix(t(a1[[i]]))
  a2[[i]] <- as.matrix(t(a2[[i]]))
  colnames(a1[[i]]) <- paste0("BL6_", proper(strains[i]),"xBl6_", colnames(a1[[i]]))
  colnames(a2[[i]]) <- paste0(strains[i], "_", proper(strains[i]),"xBl6_", colnames(a2[[i]]))
  count_mat[[i]] <- aggregate.Matrix(t(assays(allelic_count_filt[[i]])[['tot']]), 
                                groupings = olData(allelic_count_filt[[i]])$seurat_clusters, fun = "sum")
  count_mat[[i]] <- as.matrix(t(count_mat[[i]]))
  df[[i]] <- data.frame(clone = colnames(count_mat[[i]]))
  dds[[i]] <- DESeqDataSetFromMatrix(countData = count_mat[[i]], colData = df[[i]], design = ~ clone)
  dds[[i]] <- estimateSizeFactors(dds[[i]])
  sizeFactor[[i]] <- dds[[i]]$sizeFactor
  #normalising the counts
  a1_norm[[i]] <- sweep(a1[[i]], 2, sizeFactor[[i]], FUN = "/")
  a2_norm[[i]] <- sweep(a2[[i]], 2, sizeFactor[[i]], FUN = "/")
  colnames(a1_norm[[i]]) <- paste0("BL6_", proper(strain),"xBl6_", colnames(a1_norm[[i]]))
  colnames(a2_norm[[i]]) <- paste0(strain, "_", proper(strain),"xBl6_", colnames(a2_norm[[i]]))
  dat[[i]] <- cbind(a1_norm[[i]], a2_norm[[i]])
}


# Preparing data for test
#generate total allelic count and allelic ratio tables
#sum first column with fifth, second with sixth, etc.
tot_counts <- list()
genes_select <- list()
a1_counts <- list()
ai_ratio <- list()
for (i in 1:length(allelic_count_filt)){
  tot_counts[[i]] <- sapply(1:(ncol(dat[[i]])/2), function(x){rowSums(dat[[i]][,c(x,x+ncol(dat[[i]])/2)])})
  colnames(tot_counts[[i]]) <- gsub("^[^_]*_", "", colnames(dat[[i]])[1:ncol(tot_counts[[i]])])
  #select genes that have at least 10 counts across all replicates
  genes_select[[i]] <- rownames(tot_counts[[i]])[rowSums(tot_counts[[i]] >= 10) == ncol(tot_counts[[i]])]
  tot_counts[[i]] <- tot_counts[[i]][genes_select[[i]],]
  a1_counts[[i]] <- dat[[i]][genes_select[[i]], 1:ncol(tot_counts[[i]])]
  ai_ratio[[i]] <- a1_counts[[i]]/tot_counts[[i]]
  colnames(ai_ratio[[i]]) <- gsub("^[^_]*_", "", colnames(dat[[i]])[1:ncol(tot_counts[[i]])])
}


# Estimating global beta-binomial distribution parameters
genesXY <- read.table(XY_genes)
genesIMPR <- read.xlsx(imprinted_genes)
global_estim <- list()
#estimating global parameters of Beta-Binomial distribution
for (i in 1:length(allelic_count_filt)){
  global_estim[[i]] <- glob_disp(a1_counts[[i]], tot_counts[[i]], genesXY, genesIMPR)
}


# Estimating beta-binomial parameters and performing dispersion shrinkage
param_estims <- list()
param_reestim <- list()
for (i in 1:length(allelic_count_filt)){
  param_estims[[i]] <- estim_params_bulk(a1_counts[[i]], tot_counts[[i]])
  param_estims[[i]] <- param_estims[[i]][rowSums(is.na(param_estims[[i]])) != ncol(param_estims[[i]]), ]
  #paramater delta regulates how tightly dispersion is shrunk towards the mean
  #change delta=2000 for SPRETxBl6 strain
  param_reestim[[i]] <- correct_theta(param_estims[[i]], delta = 1000)
}

# Running beta-binomial test
bb_test_res <- list()
for (i in 1:length(allelic_count_filt)){
  bb_test_res[[i]] <- beta_binom_test_bulk(a1_counts[[i]], tot_counts[[i]], param_reestim[[i]], global_estim[[i]])
  #without shrinkage
  bb_test_res[[i]]$fdr_orig <- p.adjust(bb_test_res[[i]]$pval_orig, method = "fdr")
  #with shrinkage
  bb_test_res[[i]]$fdr_shrunk <- p.adjust(bb_test_res[[i]]$pval_adj, method = "fdr")
  bb_test_res[[i]] <- bb_test_res[[i]][order(bb_test_res[[i]]$fdr_shrunk),]
}

test_names <- c("CastxBl6", "MolfxBl6", "PwkxBl6", "SpretxBl6")

filename <- "/data/bbtest_res_ASPEN/bbtest_scRNA_exons_dispshrunk.xlsx"
excel <- createWorkbook(filename)
for (i in 1:length(bb_test_res)){
  addWorksheet(excel, test_names[i])
  writeData(excel, sheet = i, bb_test_res[[i]], rowNames=T, colNames=T)
  saveWorkbook(excel, file = filename, overwrite = TRUE)
}



