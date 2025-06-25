# Bradley June 2025 #

# Specify module
# 
library(Seurat)

# Specify if testing
testing=T

# Load in Seurat object and calculate IFN burden score
args = commandArgs(trailingOnly=TRUE)
seurf=args[1] # seuat file.
sampcol = args[2] # Sample column
indcol=args[3] # colname that specifies individual label.
ctcol=args[4] # colname that specifies column name.
gene_listf=args[5] # tab delim list of genes
gene_list_name = args[6]
pathout=args[7] # pathout
if(testing){
    seurf = "input/"
    sampcol = "sanger_sample_id"
    indcol = "Genotyping_ID"
    ctcol = "predicted_labels"
    gene_listf = "input/ISG_fixed_all.txt"
    gene_list_name = "ISG"
    outname = "results/"
}


# 1. Load object
print("..Loading in object")
seurat_obj <- readRDS(seurf)

# 2. Load gene set
gene_list = list(gene_list_name = read.delim(gene_listf, sep = "\t"))

# 3. Calculate score
print("..Calculating scores")
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = gene_list_ngene_listamed,
  name = gene_list_name
)

# 4. Extract and save the output
to_save = seurat_obj@meta.data[, c(sampcol, indcol, ctcol, gene_list_name)]

if(!file.exists(outname)){
    dir.create(outname)
}

fullout = paste0(outname, gene_list_name, "-module_scores.txt.gz")

print("..Saving")
write.table(
  to_save,
  file = gzfile(fullout),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)






