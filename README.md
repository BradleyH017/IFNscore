# Code used to estimate IFN module scores in IBDverse (and subsequent analysis)

### Installation

If using Sanger farm. No downloads or installs are required. 
Otherwise, singularity image (HGI/softpack/users/eh19/test-scvi-reserve/33) can be downloaded from dockerhub:
```
export SINGULARITY_CACHEDIR=$PWD/.singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
singularity pull docker://bh18/atlassing:33
```

#### 1. Clean anndata, convert to seurat
- Remove the .var info (make sure gene names are in index and are unique)
```
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
module load HGI/softpack/users/eh19/test-scvi-reserve/33 # Or run in singalarity container installed above
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
- Using the Seurat object in R
```
mkdir -p logs
MEM=500000
bsub -J "get_module_scores" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/get_module_scores-stderr \
    -o logs/get_module_scores-stdout \
    "Rscript scripts/estimate_ModuleScores.r 'input/seurat_obj.rds' 'input/gene_list/schoggins_379ISGs.txt' 'schogginsISG' 'results/'"
```

- Or using the scanpy object in python (if scaling '-s' - requires lots of memory. Otherwise ~350GB will do it)
```
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
MEM=1000000
module load HGI/softpack/users/eh19/test-scvi-reserve/33
bsub -J "get_module_scores_py" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 -q teramem \
    -e logs/get_module_scores-%J-stderr \
    -o logs/get_module_scores-%J-stdout \
    "python scripts/estimate_ModuleScores.py -i ${input_f} -m "sanger_sample_id,Genotyping_ID,disease_status,predicted_category,tissue,predicted_labels" -gl 'input/gene_list/schoggins_379ISGs_fixed.txt' -gn 'schogginsISG_fixed' -o 'results/'"
```