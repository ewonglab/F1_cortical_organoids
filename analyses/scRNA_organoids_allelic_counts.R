

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

allele_split <- function(obj){
  
  #Different versions of SingleCellExperiment
  #a1 <- round(obj@assays[["RNA"]]@data[!grepl("pseudo", rownames(obj@assays[["RNA"]]@data)),])
  #a2 <- round(obj@assays[["RNA"]]@data[grepl("pseudo", rownames(obj@assays[["RNA"]]@data)),])
  
  a1 <- round(LayerData(obj)[!grepl("pseudo", rownames(LayerData(obj))),])
  a2 <- round(LayerData(obj)[grepl("pseudo", rownames( LayerData(obj))),])
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




# scRNA libraries were multiplexed together by clone from each genetic background. 
# The libraries were demultiplexed by genotype and mapped separately using STARSolo,
# resulting in four clones across four genetic backgrounds.
# Start by importing data for each clone as a Seurat object and select cell barcodes that passed filter using non-allelic counts

strains <- c("CAST", "MOLF", "PWK", "SPRET")
replicates <- c("clone1", "clone2", "clone3", "clone4")

# Loading barcodes of the cells that passed filter using non-allelic counts
annot <- list()
for (i in 1:length(replicates)){
  annot[[i]] <- read_xlsx("/data/To_share_Good_cells_ALL.xlsx", sheet = i+1, col_names = T)
}

# Loading STARSolo output
import_seurat <- function(dir, project){
  filt.data <- Read10X(dir)
  seur.data <- CreateSeuratObject(counts = filt.data, project = project, min.cells = 5, min.features = 100)
}  


# Clone 1 
res_dir <- "/align/GRCm38/scRNA_organoids/"
filt_dir <- paste0(res_dir, "clone1/", strains, "/STARsoloIndex/clone1_", strains, "_STARsoloIndex_Solo.out/Gene/filtered")

clone1.seur.obj <- sapply(filt_dir, import_seurat, project = replicates[1])
names(clone1.seur.obj) <- paste0(strains, "xBL6")

#excluding cells that did not pass the filter
clone1.seur.sub <- list()
for (i in 1:length(replicates)){
  clone1.seur.sub[[i]] <- clone1.seur.obj[[i]][,rownames(clone1.seur.obj[[i]]@meta.data) %in% annot[[1]]$cell_barcode_clean]
}


# Clone 2
filt_dir <- paste0(res_dir, "clone2/", strains, "/STARsoloIndex/clone2_", strains, "_STARsoloIndex_Solo.out/Gene/filtered")

clone2.seur.obj <- sapply(filt_dir, import_seurat, project = replicates[2])
names(clone2.seur.obj) <- paste0(strains, "xBL6")

clone2.seur.sub <- list()
for (i in 1:length(replicates)){
  clone2.seur.sub[[i]] <- clone2.seur.obj[[i]][,rownames(clone2.seur.obj[[i]]@meta.data) %in% annot[[2]]$cell_barcode_clean]
}


# Clone 3
filt_dir <- paste0(res_dir, "clone3/", strains, "/STARsoloIndex/clone3_", strains, "_STARsoloIndex_Solo.out/Gene/filtered")

clone3.seur.obj <- sapply(filt_dir, import_seurat, project = replicates[3])
names(clone3.seur.obj) <- paste0(strains, "xBL6")

clone3.seur.sub <- list()
for (i in 1:length(replicates)){
  clone3.seur.sub[[i]] <- clone3.seur.obj[[i]][,rownames(clone3.seur.obj[[i]]@meta.data) %in% annot[[3]]$cell_barcode_clean]
}


# Clone 4
filt_dir <- paste0(res_dir, "clone4/", strains, "/STARsoloIndex/clone4_", strains, "_STARsoloIndex_Solo.out/Gene/filtered")

clone4.seur.obj <- sapply(filt_dir, import_seurat, project = replicates[4])
names(clone4.seur.obj) <- paste0(strains, "xBL6")

clone4.seur.sub <- list()
for (i in 1:length(replicates)){
  clone4.seur.sub[[i]] <- clone4.seur.obj[[i]][,rownames(clone4.seur.obj[[i]]@meta.data) %in% annot[[4]]$cell_barcode_clean]
}


# Data for individual clones is merged by strain
CAST.seur.mrg <- merge(clone1.seur.sub[[1]], y = c(clone2.seur.sub[[1]], 
                                                   clone3.seur.sub[[1]], 
                                                   clone4.seur.sub[[1]]), 
                       add.cell.ids = replicates, project = "CASTxBl6")

MOLF.seur.mrg <- merge(clone1.seur.sub[[2]], y = c(clone2.seur.sub[[2]], 
                                                   clone3.seur.sub[[2]], 
                                                   clone4.seur.sub[[2]]), 
                       add.cell.ids = replicates, project = "MOLFxBl6")

PWK.seur.mrg <- merge(clone1.seur.sub[[3]], y = c(clone2.seur.sub[[3]], 
                                                  clone3.seur.sub[[3]], 
                                                  clone4.seur.sub[[3]]), 
                      add.cell.ids = replicates, project = "PWKxBl6")

SPRET.seur.mrg <- merge(clone1.seur.sub[[4]], y = c(clone2.seur.sub[[4]], 
                                                    clone3.seur.sub[[4]], 
                                                    clone4.seur.sub[[4]]), 
                        add.cell.ids = replicates, project = "SPRETxBl6")


# Splitting STARSolo output by alleles
CAST.organoids.allelic.counts <- allele_split(CAST.seur.mrg)
MOLF.organoids.allelic.counts <- allele_split(MOLF.seur.mrg)
PWK.organoids.allelic.counts <- allele_split(PWK.seur.mrg)
SPRET.organoids.allelic.counts <- allele_split(SPRET.seur.mrg)

saveRDS(CAST.organoids.allelic.counts, paste0(dir, "CAST.organoids.allelic.counts.exons.rds"))
saveRDS(MOLF.organoids.allelic.counts, paste0(dir, "MOLF.organoids.allelic.counts.exons.rds"))
saveRDS(PWK.organoids.allelic.counts, paste0(dir, "PWK.organoids.allelic.counts.exons.rds"))
saveRDS(SPRET.organoids.allelic.counts, paste0(dir, "SPRET.organoids.allelic.counts.exons.rds"))





