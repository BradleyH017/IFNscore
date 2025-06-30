# Bradley June 2025 #

# Specify module
library(Seurat)
library(SeuratDisk)
library(tidyverse)

# Load in Seurat object and calculate IFN burden score
args = commandArgs(trailingOnly=TRUE)
seurf=args[1] # seurf = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/edQTLs/IFNscore/input/subsampled.h5seurat"
gene_listf=args[3] # "input/schoggins_379ISGs.txt"
gene_list_name = args[4] # gene_list_name="ISG"
outname=args[5] # outname="results/"
metadata_cols = unlist(strsplit(metadata_cols, ","))

# 1. Load object
print("..Loading in object")
seurat_obj <- readRDS(seurf)

# 2. Load gene set
gene_vector <- read.delim(gene_listf, sep = "\t", header = F)[,1]
gene_list <- list()
gene_list[[gene_list_name]] <- gene_vector

# 3. Calculate score
print("..Calculating scores")
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = gene_list,
  name = gene_list_name
)

# Rename (suffix added)
old = paste0(gene_list_name, "1")
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == old)] = gene_list_name

# 4. Extract and save the output
to_save = seurat_obj@meta.data %>% 
  rownames_to_column(var = "cell") %>% 
  select(cell, !!gene_list_name)

if(!file.exists(outname)){
    dir.create(outname)
}

fullout = paste0(outname, gene_list_name, "-module_scores_seurat.txt.gz")

print("..Saving")
write.table(
  to_save,
  file = gzfile(fullout),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)






