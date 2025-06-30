# Code used to estimate IFN module scores in IBDverse (and subsequent analysis)

#### 1: Make seurat from expression and metadata manually
Save a list of samples from the anndata object. The -l option means the matrix being saved in the log1p_cp10k matrix
```
module load HGI/softpack/users/eh19/test-scvi-reserve/33
MEM=50000
input_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
bsub -J "get_samp_list_chunk" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/get_samp_list_chunk-%J-stderr \
    -o logs/get_samp_list_chunk-%J-stdout \
    "python scripts/get_samp_list_chunk.py -i ${input_f} -l 'log1p_cp10k' -s 'sanger_sample_id' -m 'sanger_sample_id,Genotyping_ID,predicted_labels,predicted_category,disease_status' -o 'input/expr_per_sample'"
```

Then run a snakemake pipeline to: 1) Extract mtx for each, and 2) Collate the mtx into a single seurat object
```
bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_CONV" < submit_snakemake_BH.sh 
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