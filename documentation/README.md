# PreciCE
## A unified workflow for data-driven precision cell fate engineering via highly multiplexed gene control


### NOTES ON THIS VERSION ### 
We have optimized PreciCE by enabling the entire workflow (including gene regulatory network inference) to be run in a single Python session, and performed additional streamlining and bug fixes, as specified in the commit comments. 


### INSTALLATION ###
**Conda (exact, reproducible environment -- Linux x86_64 only for now)**

This integrated workflow is sensitive to dependency versions. We  recommend recreating the environment from the explicit Conda spec:

```bash
# from a shell in /path/to/precice/documentation (Linux x86_64)
conda create -n environment_name --file conda-spec-linux-64.txt
conda activate environment_name
```

**Downloading large input files:**
Download the following (human-specific) files to /path/to/precice/input/resources. Note: use these v9 files, not the newer v10 files unless you specifically plan to work with SCENIC+ rather than SCENIC, as specified by the Aerts lab.

hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather (https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather)

motifs-v9-nr.hgnc-m0.001-o0.0.tbl (https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl)

Optional: Download the example scRNA-seq dataset from Friedman et al., Cell Stem Cell (2018), pre-formatted as an .h5ad file for use with PreciCE. This dataset is provided as part of the GitHub release assets.


### GENERAL NOTES ###
## Now uses different data processing steps for differential expression and GRN inference
The workflow now performs highly-variable-gene filtering for differential expression but no such filtering for GRN inference. The pipeline will currently output two separate processed datasets - one with HVG filtering and the other without. The pipeline will import the relevant version before running DE and GRN inference, respectively, by redefining "adata" right before each function is run. Make sure the right version of adata is used (e.g., ending in ...for_DE.csv for differential expression; ...for_pyscenic.csv for everything else.)

## Recommended parameter sweep setting for GRN inference
Gene regulatory network inference is sensitive to parameters in pySCENIC's prune2df() function (in pyscenic_pipeline.py). We recommend running a parameter sweep, i.e., rerunning the pipeline with different settings for prune2df() to gauge the robustness of the perturbation in the final output plot.
We recommend the following settings for a prune2df() parameter sweep

# prune2df() permissiveness setting 1 (pySCENIC python default):
rank_threshold=1500, auc_threshold=0.05, nes_threshold=3.0, motif_similarity_fdr=0.001

# prune2df() permissiveness setting 2 (More permissive; pySCENIC command line interface default):
rank_threshold=5000, auc_threshold=0.05, nes_threshold=3.0, motif_similarity_fdr=0.001

# prune2df() permissiveness setting 3 (Very permissive):
rank_threshold=5000, auc_threshold=0.03, nes_threshold=2.0, motif_similarity_fdr=0.01

# prune2df() permissiveness setting 4 (Recklessly permissive -- will likely include many false-positive network edges):
rank_threshold=5000, auc_threshold=0.01, nes_threshold=1.0, motif_similarity_fdr=0.01


## Increasing resolution in the final output plot
When running run_precice(), sometimes two genes are ranked equally highly and show up as two genes per row in the final output plot. To resolve these, change these parameters lambda_max_step and/or lambda_step_size in run_precice (Note: significantly increases compute time).


### LAST UPDATE: JENS MAGNUSSON 2025-09-30 ###