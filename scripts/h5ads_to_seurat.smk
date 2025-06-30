# Get the gene list
with open(config["sample_list"], "r") as f:
    samples = [line.strip() for line in f if line.strip()]

# For each gene, run the clumping
rule run_all:
	input:
		expand("{outdir}/../seurat_all.rds", outdir=config["outdir"])
    
rule h5d_to_mtx:
    output:
        "{outdir}/{sample}/matrix.mtx.gz" 
    params:
        outdir=config["outdir"]
    threads: 1
    resources: 
        mem_mb=8000,
        disk_mb=8000,
        queue="normal"
    singularity:
        "/software/hgi/softpack/installs/users/eh19//test-scvi-reserve/33-scripts/singularity.sif"
    shell:
        r"""
        python scripts/h5_to_mtx.py '{params.outdir}/{wildcards.sample}/{wildcards.sample}.h5ad' '{params.outdir}/{wildcards.sample}/' --out_file '' --verbose
        """

def gather_per_sample_outputs(wildcards):
    return expand("{outdir}/{sample}/matrix.mtx.gz",
                  outdir=config["outdir"], sample=samples)

rule aggregate:
    input:
        gather_per_sample_outputs
    output:
        "{outdir}/../seurat_all.rds"
    params:
        outdir=config["outdir"],
        sample_list=config["sample_list"]
    threads: 1
    resources:
        mem_mb=750000,
        disk_mb=750000,
        queue='teramem',
        threads=1
    singularity:
        "/software/hgi/softpack/installs/users/tr12//JAGUAR_CytoOmic/2-scripts/singularity.sif"
    shell:
        r"""
        Rscript scripts/mtx_to_seurat.r {params.sample_list} '{params.outdir}/' '{params.outdir}/../'
        """