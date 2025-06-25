# Code used to estimate IFN module scores in IBDverse (and subsequent analysis)

#### 1. Convert the h5ad object to seurat
```
mkdir -p logs
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
output_f="input/celltypist_0.5_ngene_ncount_mt_filt_nomiss."
module load HGI/common/h5ad_to_seurat/0.1
MEM=150000
bsub -J "conv_seurat" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/conv_seurat-stderr \
    -o logs/conv_seurat-stdout \
    "h5ad_to_seurat ${input_f} ${output_f}"
```

Also add the module genes to each line of tab delim file. Put this in "input/"

#### 2. Calculate the module scores for each cell-type, within each sample
```
mkdir -p logs
MEM=500000
bsub -J "get_module_scores" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/get_module_scores-stderr \
    -o logs/get_module_scores-stdout \
    "Rscript scripts/estimate_ModuleScores.r 'input/' 'sanger_sample_id' 'Genotyping_ID' 'predicted_labels' 'input/ISG_fixed_all.txt' 'ISG' 'results/'    > \
    logs/get_module_scores-${level}-${hm}.Rout"
```