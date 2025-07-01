# Code used to estimate IFN module scores in IBDverse (and subsequent analysis)

#### 1. Clean anndata, convert to seurat
- Remove the .var info (make sure gene names are in index and are unique)
```
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
module load HGI/softpack/users/eh19/test-scvi-reserve/33
python scripts/clean_anndata.py -i ${input_f} -l 'log1p_cp10k' -m "sanger_sample_id,Genotyping_ID,predicted_labels,predicted_category,disease_status" -o "input/"
```

- Now submit the conversion job
```
MEM=500000
module load HGI/common/h5ad_to_seurat/0.1
mkdir -p logs
bsub -J "conv_seurat" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
        -e logs/conv_seurat-%J-stderr \
        -o logs/conv_seurat-%J-stdout \
    "h5ad_to_seurat input/clean.h5ad"

mv clean.h5seurat input/
```

#### 2. Calculate the module scores for each cell-type, within each sample
```
mkdir -p logs
MEM=500000
bsub -J "get_module_scores" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/get_module_scores-stderr \
    -o logs/get_module_scores-stdout \
    "Rscript scripts/estimate_ModuleScores.r 'input/seurat_obj.rds' 'input/gene_list/schoggins_379ISGs.txt' 'schogginsISG' 'results/'"
```