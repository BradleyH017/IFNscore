# Code used to estimate IFN module scores in IBDverse (and subsequent analysis)

#### 1. Convert the h5ad object to seurat
```
module load HGI/common/h5ad_to_seurat/0.1
mkdir -p input
h5ad_to_seurat /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad input/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad
```