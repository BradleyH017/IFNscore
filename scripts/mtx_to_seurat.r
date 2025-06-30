# Bradley June 2025
# Creating seurat object from mtx 
# module load HGI/softpack/users/tr12/JAGUAR_CytoOmic/2

library(Seurat)

args = commandArgs(trailingOnly=TRUE)
sample_list <- args[1] # sample_list="input/expr_per_sample/samples.txt"
input_dir <- args[2] # input_dir="input/expr_per_sample/"
out_dir <- args[3] # out_dir="input/"

# Load sample_list
print("..Loading in per sample matrices")
samples = readLines(sample_list)
samp_list = vector("list", length = length(samples))
for(s in samples){
    print(s)
    expression_data <- Read10X(data.dir = paste0(input_dir, s))
    seurat_obj <- CreateSeuratObject(counts = expression_data, project = s)
    seurat_obj@assays$RNA$data <- expression_data
    samp_list[[which(samples == s)]] <- seurat_obj
}

# Combine
print("Combining")
merged_seurat <- merge(
  x = samp_list[[1]], 
  y = samp_list[2:length(samp_list)]
)

# Joing layers
print("Merging layers")
merged_seurat <- JoinLayers(merged_seurat)

# Save
print("..Saving")
saveRDS(seurat_obj, file = paste0(out_dir, "seurat_obj.rds"))
