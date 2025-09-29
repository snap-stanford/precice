
### RECOMMENDED WORKFLOW FOR RUNNING PRECICE IN PYTHON ###
# Includes data preprocessing, differential expression, gene regulatory network inference, gene prediction #

# In Python v3.10.17:
import os
os.chdir('/path/to/precice/model')  # Set default working directory (currently required for correct paths)

import precice                      # PreciCE main functions
import create_network               # For visualizing network with third-party tools
import learn_weights                # For tuning network weights
import plot                         # For plotting results
import precice_utils                # Helper functions
import pyscenic_utils               # NEW! Runs network inference in Python session
import run_precice                  # PreciCE helper functions
import utils                        # Helper functions
import scanpy as sc
import sys
from plot import make_precision_plot    # For plotting precice output
from types import SimpleNamespace       # For SimpleNamespace()
from pathlib import Path                # Required for creating and using project_directory

### 1. DEFINE KEY VARIABLES ###
project_directory = Path("../data/output/2025-09-29_FriedmanStemToMeso")
name = "Friedman.h5ad"          # Must match the .h5ad file name
input_path = '../data/input/datasets/Friedman.h5ad'
raw_adata = sc.read_h5ad('../data/input/datasets/Friedman.h5ad')

source_name = 'stem'
target_name = 'meso'
cell_type_label = 'label'
batch_key = 'day'

### 2. INITIALIZE PRECICE ###
workflow = precice.precice(adata=raw_adata, 
                   path = input_path, 
                   dir = project_directory,
                   cell_filter=True,
                   gene_filter=True,
                   mito_percent_max=5) # Set permissiveness for mitochondrial read inclusion (e.g., default: 5, permissive: 12)

### 3. DIFFERENTIAL EXPRESSION ###
# Read the adata file that was generated for differential expression in step 2 (WITH filter for highly variable genes!)
de_path = project_directory / f"{name}_processed_for_DE.h5ad"
workflow.adata = sc.read_h5ad(de_path)
# Print info statement specifying what is being passed as adata (sanity check to verify that the correct data version is being passed):
print(                                             # INFO message
    f"[INFO] Differential expression being run on "
    f"{de_path.name} [{workflow.adata.n_obs} cells × {workflow.adata.n_vars} genes]"
)
# Perform differential expression and save file to project_directory:
workflow.set_up_scvi(batch_key=batch_key)
de_table = workflow.get_DE(
    source_name=source_name,                  # Starting cell type
    target_name=target_name,                  # Target cell type
    cell_type_label=cell_type_label,
    batch_key=batch_key,
)
workflow.write_DE_files(DE_dir=project_directory) # Export to file

### 4. GENE REGULATORY NETWORK INFERENCE ###
# Read the adata file that was previously generated for pySCENIC in step 2 (NOT FILTERED for highly variable genes!)
ps_path = project_directory / f"{name}_processed_for_pyscenic.h5ad"
workflow.adata = sc.read_h5ad(ps_path)
# Print info statement specifying what is being passed as adata (sanity check to verify that the correct data version is being passed)
print(                                             # INFO message
    f"[INFO] pySCENIC being run on "
    f"{ps_path.name} [{workflow.adata.n_obs} cells × {workflow.adata.n_vars} genes]"
)
# Run pyscenic to compute GRN
workflow.run_pyscenic(
    species      = "human",            # e.g. 'human' or 'mouse'
    metadata_key = cell_type_label,    # your obs column for cell‐type subsetting
    n_workers    = 36,                 # 64 Dask workers ⇢ 16 CPU cores
    seed         = 777,                # for reproducibility
    out_dir      = project_directory,  
)
# Learn edge weights:
weights_args = SimpleNamespace(
    data_path  = project_directory / f"{name}_processed_for_pyscenic.h5ad",
    graph_name = project_directory / f"{name}_allcells_linearweights.csv",
    out_dir    = project_directory,
    method     = "linear"
)
learn_weights.main(weights_args)


# Import network (either generated above or pre-existing) for downstream use. Bypasses old code that relied on an intermediate .csv file that stored file paths but is back-compatible with that option:
network_file = project_directory / "Friedman_allcells_linearweights_linear_learntweights.csv" #customize this to specify the actual network file
workflow.get_network(network_path=network_file)
transition = source_name +'_to_' + target_name
python_path = sys.executable

### 5. RUN PRECICE TO COMPUTE PERTURBATIONS ###
workflow.run_precice(
    species       = "human",
    python_path   = sys.executable,
    network_path  = workflow.network_path,
    DE_path       = workflow.DE_filenames[transition],
    out_dir       = project_directory,
    lambda_max_step  = 1,       # ← Default 1. Increase to resolve two-genes-per row in plot's first row (e.g., change to 10)
    lambda_step_size = 1e-4,     # ← Default 1e-4. Decrease to resolve two-genes-per row in plot's subsequent rows (e.g. change to 1e-5)
)

### 6. PLOT RESULTS ###
plot_input_file = project_directory / "active_no_parallel_p_0.02_max_step_1_hops_3_DE_DE_Friedman_stem_to_meso_adjacency_Friedman_allcells_linearweights_linear_learntweights.pkl" #.pkl file generated by run_precice()
plot_output_file = "2025-09-29_FriedmanStemToMeso.pdf" # Optional: specify .jpg file type

make_precision_plot(
    fname        = plot_input_file,
    max_rows     = 12,
    source       = "Stem",
    target       = "Meso",
    out_dir      = project_directory,
    out_name     = plot_output_file,
)

### 7. (OPTIONAL) EXPORT GRN FOR VISUALIZATION WITH EXTERNAL TOOLS. NOTE: NOT FUNCTIONAL YET (2025-09-29) ###
#Export network for visualization (e.g., with Observable-hq)
regulon_path = project_directory / "allcells"   #Path to regulons folder (where each file is one regulon)
args = SimpleNamespace(
    regulons = str(regulon_path),  # ← convert Path → str (or use regulon_path.as_posix())
    species  = "human",
    name     = None                # let create_net() derive a name from the folder
)
create_network.create_net(args)





