# Bradley June 2025
# Creating seurat object from mtx 
# module load HGI/softpack/users/tr12/JAGUAR_CytoOmic/2

library(Seurat)

args = commandArgs(trailingOnly=TRUE)
data_dir <- args[1] # data_dir="input/mtx"
out_dir <- args[2] # input/

# Load expression
print("..Reading 10X")
expression_data <- Read10X(data.dir = data_dir)

# Create Seurat object
print("Creating seurat object")
seurat_obj <- CreateSeuratObject(counts = expression_data, project = "MyProject")
seurat_obj@assays$RNA$data <- expression_data


# Save
print("..Saving")
saveRDS(seurat_obj, file = paste0(out_dir, "seurat_obj.rds"))
