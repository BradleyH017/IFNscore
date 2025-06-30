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
            '-s', '--sample_col',
            action='store',
            dest='metadata_cols',
            required=True,
            help=''
        )

    parser.add_argument(
            '-l', '--layer',
            action='store',
            dest='layer',
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
            '-o', '--outdir',
            action='store',
            dest='outdir',
            required=True,
            help=''
        )
    
    return parser.parse_args()


def main():
    # Parse options
    inherited_options = parse_options()
    
    # Load object
    print("..Loading h5ad")
    adata = sc.read_h5ad(inherited_options.input_file)
    
    # Make sure .X is the correct layer
    adata.X = adata.layers[inherited_options.layer]
    
    # Make sure outdir exists
    os.makedirs(inherited_options.outdir, exist_ok=True)
    
    # Save a list of unique values for chunks
    chunks = adata.obs[inherited_options.sample_col].reset_index(drop=True).drop_duplicates().tolist()
    with open(f"{inherited_options.outdir}/samples.txt", "w") as f:
        for item in chunks:
            f.write(f"{item}\n")
            
    # Save a .h5ad per sample (in new dir)
    for item in chunks:
        print(item)
        outdir=f"{inherited_options.outdir}/{item}"    
        os.makedirs(outdir, exist_ok=True)
        adata[adata.obs[inherited_options.sample_col] == item].write_h5ad(f"{outdir}/{item}.h5ad")
        
    # Save the metadata
    to_save = adata.obs[inherited_options.metadata_cols.split(",")]   
    to_save.to_csv(f"{inherited_options.outdir}/metadata.txt.gz", sep="\t", index=True, compression='gzip')
        
    
    
