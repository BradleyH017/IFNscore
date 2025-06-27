# Code used to estimate IFN module scores in IBDverse (and subsequent analysis)

#### 1. Convert the h5ad object to seurat
```
mkdir -p logs
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
output_f="input/celltypist_0.5_ngene_ncount_mt_filt_nomiss."
module load HGI/common/h5ad_to_seurat/0.1
MEM=50000
bsub -J "conv_seurat" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/conv_seurat-%J-stderr \
    -o logs/conv_seurat-%J-stdout \
    "h5ad_to_seurat ${input_f} ${output_f}"
```

#### 1-Alternative: Make seurat from expression and metadata manually
Extract this from the h5ad and save as mtx
```
module load HGI/softpack/users/eh19/test-scvi-reserve/33
mkdir -p logs
mkdir -p input/mtx
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
MEM=500000
bsub -J "h5_to_mtx" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/h5_to_mtx-%J-stderr \
    -o logs/h5_to_mtx-%J-stdout \
    "python scripts/h5_to_mtx.py ${input_f} 'input/mtx' --out_file '' --verbose"
```

Now make a seurat object from this
```
module load HGI/softpack/users/tr12/JAGUAR_CytoOmic/2
MEM=500000
bsub -J "mtx_to_seurat" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/mtx_to_seurat-%J-stderr \
    -o logs/mtx_to_seurat-%J-stdout \
    "Rscript scripts/mtx_to_seurat.py 'input/mtx' 'input/'"
```


#### 2. Calculate the module scores for each cell-type, within each sample
```
mkdir -p logs
MEM=500000
bsub -J "get_module_scores" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/get_module_scores-stderr \
    -o logs/get_module_scores-stdout \
    "Rscript scripts/estimate_ModuleScores.r 'input/seurat_obj.rds' 'input/schoggins_379ISGs.txt' 'schogginsISG' 'results/'"
```