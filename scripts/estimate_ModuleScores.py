
__author__ = 'Bradley Harris'
__date__ = '2025-06-25'
__version__ = '0.0.1'

# Change dir
import os
cwd = os.getcwd()
print(cwd)

# Load packages
import sys
import os
sys.path.append('/software/team152/bh18/pip')
sys.path.append('/usr/local/')
print("System path")
print(sys.path)
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mp
from matplotlib import pyplot as plt
from matplotlib.pyplot import rc_context

def parse_options():    
    # Inherit options
    parser = argparse.ArgumentParser(
            description="""
                Estimating module scores for a specific set gene set
                """
        )
    
    parser.add_argument(
            '-i', '--input_file',
            action='store',
            dest='input_file',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-m', '--metadata_cols',
            action='store',
            dest='metadata_cols',
            required=True,
            help=''
        )

    parser.add_argument(
            '-gl', '--gene_listf',
            action='store',
            dest='gene_listf',
            required=True,
            help=''
        )

    parser.add_argument(
            '-gn', '--gene_list_name',
            action='store',
            dest='gene_list_name',
            required=True,
            help=''
        )

    parser.add_argument(
            '-o', '--outname',
            action='store',
            dest='outname',
            required=True,
            help=''
        )

    return parser.parse_args()


def main():
    # Parse options
    inherited_options = parse_options()
    input_file = inherited_options.input_file # input_file="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_0.5_ngene_ncount_mt_filt_nomiss.h5ad"
    metadata_cols = inherited_options.metadata_cols # metadata_cols = "sanger_sample_id,Genotyping_ID,disease_status,predicted_category,tissue,predicted_labels"
    metadata_cols = metadata_cols.split(",")
    gene_listf = inherited_options.gene_listf # gene_listf = "input/schoggins_379ISGs.txt"
    gene_list_name = inherited_options.gene_list_name # gene_list_name = "379ISGs"
    outname = inherited_options.outname # outname = "results/"
    
    # 1. Load
    print("..Loading object")
    adata = sc.read_h5ad(input_file)
    
    # 2. Load gene set
    with open(gene_listf) as f:
        gene_list = [line.strip() for line in f if line.strip()]
        
    # 3. Compute module scores using scanpy
    # Replace ENS with gene names if these are used as rows
    if adata.var_names[0].startswith("ENS"):
        adata.var['ENS'] = adata.var_names.copy()
        adata.var_names = list(adata.var['gene_symbols'])
        adata.var_names_make_unique()
        adata.var['gene_symbols_unique'] = adata.var_names.copy()
        
    # make sure .X is normalised data
    adata.X = adata.layers['log1p_cp10k']
    
    print("..Computing module scores")
    sc.tl.score_genes(adata, gene_list, score_name=gene_list_name)
    
    # 4. Extract and save
    print("..Saving")
    metadata_cols.append(gene_list_name)
    to_save = adata.obs[metadata_cols]   
    os.makedirs(outname, exist_ok=True)
    fullout = f"{outname}{gene_list_name}-module_scores.txt.gz"
    to_save.to_csv(fullout, sep="\t", index=True, compression='gzip')
    
    
# Execute
if __name__ == '__main__':
    main()

