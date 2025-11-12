# PreciCE
## A unified workflow for data-driven precision cell fate engineering via highly multiplexed gene control. 

## PreciCE will run:
1. **Data preprocessing** - Filters cells/genes, normalization
2. **Differential expression** - Identifies DE genes between cell types
3. **Network inference** - Constructs GRN using pySCENIC
4. **Perturbation analysis** - Computes gene perturbations
5. **Visualization** - Generates precision plots

### NOTES ON THIS VERSION ### 
We have optimized PreciCE by enabling the entire workflow (including gene regulatory network inference) to be run in a single Python session, and performed additional streamlining and bug fixes, as specified in the GitHub commit comments. 


### INSTALLATION ###
**Conda (exact, reproducible environment -- Linux x86_64 only for now)**

This integrated workflow is sensitive to dependency versions. We  recommend recreating the environment from the explicit Conda spec (see precice/documentation folder):

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


### USAGE ###
The PreciCE algorithm can now be run either in script mode or in interactive mode (Python REPL, for line-by-line troubleshooting), as specified below.

# Running PreciCE in script mode
### 1. Prerequisites
- Python 3.10.17
- PreciCE installed with all dependencies (see INSTALLATION above)
- Input: `.h5ad` file with cell type and batch metadata; input files specified in INSTALLATION (above)

### 2. Basic Usage
```bash
python precice_script.py \
    --work-dir /path/to/precice/model \
    --input-path /path/to/Friedman.h5ad \
    --project-dir /path/to/output/2025-01-01_FriedmanStemToMeso \
    --source-name stem \
    --target-name meso \
    --plot_output FriedmanStemToMeso.pdf
```

### 3. Input Options
## Common Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--source-name` | `stem` | Starting cell type |
| `--target-name` | `meso` | Target cell type |
| `--plot_output` | `PrecisionPlot` | File name for precision plot |
| `--cell-type-label` | `label` | Cell type column in adata.obs |
| `--batch-key` | `day` | Batch column in adata.obs |
| `--species` | `human` | Species (`human` or `mouse`) |
| `--n-workers` | `36` | Number of parallel workers |
| `--mito-percent-max` | `5` | Max mitochondrial % (use 12 for permissive) |

## Examples
**Basic run:**
```bash
python precice_script.py \
    --input-path Friedmn.h5ad \
    --source-name stem \
    --target-name meso
```

**Fine-tune plot resolution:**
```bash
python precice_script.py \
    --input-path Friedman.h5ad \
    --lambda-max-step 10 \
    --lambda-step-size 1e-5 \
    --max-rows 20
```

## Output Files
Your `--project-dir` will contain:
- `*_processed_for_DE.h5ad` - Preprocessed data for differential expression (highly-variable-gene filtering)
- `*_processed_for_pyscenic.h5ad` - Preprocessed data for pyscenic (no highly-variable-gene filtering)
- `*.pkl` - Inferred and processed network
- `allcells` - Folder containing individual regulons for selected cells (e.g., allcells)
- `DE_*_source_to_target.csv` - Differentially expressed genes
- `*_linear_learntweights.csv` - Inferred gene regulatory network
- `*_linearweights.csv` - Intermediate network file
- `*_regulons.gmt` - Regulons
- `perturbation_matrix*` - Perturbtations as matrix
- `results.pdf` - **Precision plot** (main output)

The precision plot shows genes ranked by importance for the cell state transition.

## Troubleshooting

**Import errors:** Ensure you're in the correct `--work-dir`

**File not found:** Use absolute paths for `--input-path` and `--project-dir`

**Memory issues:** Reduce `--n-workers` or use machine with more RAM

**Wrong column names:** Check your data and set `--cell-type-label` and `--batch-key`:
```python
import scanpy as sc
adata = sc.read_h5ad("data.h5ad")
print(adata.obs.columns)
```

**Runtime:** Workflow takes 2-6 hours depending on dataset size


# Running PreciCE in interactive mode (REPL):
Launch Python in your terminal:
```
python
```
Open PreciCE_Workflow_REPL.md and execute its commands line by line at the >>> prompt.


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


### LAST UPDATE: JENS MAGNUSSON 2025-11-12 ###