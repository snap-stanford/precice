#!/usr/bin/env python3
"""
PreciCE Complete Workflow Script
Includes data preprocessing, differential expression, gene regulatory network inference, and gene prediction
Python 3.10.17

Usage:
    python precice_workflow.py [options]

Examples:
    # Use all defaults
    python precice_workflow.py
    
    # Specify custom paths
    python precice_workflow.py --work-dir /path/to/precice/model --input-path /data/my_data.h5ad
    
    # Override specific parameters
    python precice_workflow.py --source-name pluripotent --target-name cardiac --n-workers 64

"""

import os
import sys
import argparse
from pathlib import Path
from types import SimpleNamespace

import scanpy as sc

# Import PreciCE modules (after setting working directory)
import precice
import create_network
import learn_weights
import plot
import precice_utils
import pyscenic_utils
import run_precice
import utils
from plot import make_precision_plot


def parse_arguments():
    """Parse command-line arguments with sensible defaults."""
    parser = argparse.ArgumentParser(
        description='Run complete PreciCE workflow for gene regulatory network analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Path arguments
    parser.add_argument('--work-dir', type=str, default='/path/to/precice/model',
                        help='Working directory for PreciCE model')
    parser.add_argument('--project-dir', type=str, default='../data/output/2025-09-24_FriedmanStemToMeso',
                        help='Project output directory')
    parser.add_argument('--input-path', type=str, default='../data/input/datasets/Friedman.h5ad',
                        help='Path to input .h5ad file')
    parser.add_argument('--name', type=str, default='Friedman.h5ad',
                        help='Name of the .h5ad file (must match input filename)')
    
    # Cell type and metadata arguments
    parser.add_argument('--source-name', type=str, default='stem',
                        help='Starting cell type name')
    parser.add_argument('--target-name', type=str, default='meso',
                        help='Target cell type name')
    parser.add_argument('--cell-type-label', type=str, default='label',
                        help='Column name for cell type labels in adata.obs')
    parser.add_argument('--batch-key', type=str, default='day',
                        help='Column name for batch information in adata.obs')
    
    # Processing parameters
    parser.add_argument('--species', type=str, default='human', choices=['human', 'mouse'],
                        help='Species for pySCENIC analysis')
    parser.add_argument('--mito-percent-max', type=float, default=5,
                        help='Maximum mitochondrial read percentage (permissiveness)')
    parser.add_argument('--n-workers', type=int, default=36,
                        help='Number of workers for parallel processing')
    parser.add_argument('--seed', type=int, default=777,
                        help='Random seed for reproducibility')
    
    # PreciCE parameters
    parser.add_argument('--lambda-max-step', type=int, default=1,
                        help='Lambda max step parameter (increase to resolve gene spacing in first row)')
    parser.add_argument('--lambda-step-size', type=float, default=1e-4,
                        help='Lambda step size parameter (decrease to resolve gene spacing in later rows)')
    
    # Plotting parameters
    parser.add_argument('--max-rows', type=int, default=12,
                        help='Maximum rows in precision plot')
    parser.add_argument('--plot-output', type=str, default='PrecisionPlot',
                        help='Output filename for plot (supports .pdf, .jpg, .png)')
    
    # Network file (optional override)
    parser.add_argument('--network-file', type=str, default=None,
                        help='Path to pre-computed network file (optional)')
    
    return parser.parse_args()


def main():
    """Main workflow for running PreciCE analysis."""
    
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set working directory (required for correct paths)
    print(f"[INFO] Setting working directory to: {args.work_dir}")
    os.chdir(args.work_dir)
    
    # =========================================================================
    # 1. DEFINE KEY VARIABLES FROM ARGUMENTS
    # =========================================================================
    project_directory = Path(args.project_dir)
    name = args.name
    input_path = args.input_path
    

    # Normalize dataset base name, so passing "Friedman.h5ad" yields "Friedman"
    # and derived files are named like "Friedman_processed_for_*.h5ad"
    dataset_stem = Path(name).stem
    
    source_name = args.source_name
    target_name = args.target_name
    cell_type_label = args.cell_type_label
    batch_key = args.batch_key
    
    print(f"[INFO] Configuration:")
    print(f"  - Project directory: {project_directory}")
    print(f"  - Input path: {input_path}")
    print(f"  - Source cell type: {source_name}")
    print(f"  - Target cell type: {target_name}")
    print(f"  - Species: {args.species}")
    
    # Read raw data
    print("[INFO] Reading raw data...")
    raw_adata = sc.read_h5ad(input_path)
    
    # =========================================================================
    # 2. INITIALIZE PRECICE
    # =========================================================================
    print("[INFO] Initializing PreciCE workflow...")
    workflow = precice.precice(
        adata=raw_adata,
        path=input_path,
        dir=project_directory,
        cell_filter=True,
        gene_filter=True,
        mito_percent_max=args.mito_percent_max
    )
    
    # =========================================================================
    # 3. DIFFERENTIAL EXPRESSION
    # =========================================================================
    print("[INFO] Running differential expression analysis...")
    
    # Read the adata file generated for DE (WITH filter for highly variable genes)
    de_path = project_directory / f"{dataset_stem}_processed_for_DE.h5ad"
    workflow.adata = sc.read_h5ad(de_path)
    
    print(
        f"[INFO] Differential expression being run on "
        f"{de_path.name} [{workflow.adata.n_obs} cells × {workflow.adata.n_vars} genes]"
    )
    
    # Perform differential expression
    workflow.set_up_scvi(batch_key=batch_key)
    de_table = workflow.get_DE(
        source_name=source_name,
        target_name=target_name,
        cell_type_label=cell_type_label,
        batch_key=batch_key,
    )
    workflow.write_DE_files(DE_dir=project_directory)
    
    # =========================================================================
    # 4. GENE REGULATORY NETWORK INFERENCE
    # =========================================================================
    print("[INFO] Running gene regulatory network inference...")
    
    # Read the adata file for pySCENIC (NOT FILTERED for highly variable genes)
    ps_path = project_directory / f"{dataset_stem}_processed_for_pyscenic.h5ad"
    workflow.adata = sc.read_h5ad(ps_path)
    
    print(
        f"[INFO] pySCENIC being run on "
        f"{ps_path.name} [{workflow.adata.n_obs} cells × {workflow.adata.n_vars} genes]"
    )
    
    # Run pySCENIC to compute GRN
    workflow.run_pyscenic(
        species=args.species,
        metadata_key=cell_type_label,
        n_workers=args.n_workers,
        seed=args.seed,
        out_dir=project_directory,
    )
    
    # Learn edge weights
    print("[INFO] Learning network edge weights...")
    weights_args = SimpleNamespace(
        data_path=project_directory / f"{dataset_stem}_processed_for_pyscenic.h5ad",
        graph_name=project_directory / f"{dataset_stem}_allcells_linearweights.csv",
        out_dir=project_directory,
        method="linear"
    )
    learn_weights.main(weights_args)
    
    # Import network for downstream use
    if args.network_file:
        network_file = Path(args.network_file)
        print(f"[INFO] Using provided network file: {network_file}")
    else:
        network_file = project_directory / f"{dataset_stem}_allcells_linearweights_linear_learntweights.csv"
        print(f"[INFO] Using default network file: {network_file}")
    
    workflow.get_network(network_path=network_file)
    transition = f"{source_name}_to_{target_name}"
    
    # =========================================================================
    # 5. RUN PRECICE TO COMPUTE PERTURBATIONS
    # =========================================================================
    print("[INFO] Running PreciCE to compute perturbations...")
    workflow.run_precice(
        species=args.species,
        python_path=sys.executable,
        network_path=workflow.network_path,
        DE_path=workflow.DE_filenames[transition],
        out_dir=project_directory,
        lambda_max_step=args.lambda_max_step,
        lambda_step_size=args.lambda_step_size,
    )
    
    # =========================================================================
    # 6. PLOT RESULTS
    # =========================================================================
    print("[INFO] Generating precision plot...")
    
    # Construct the expected pickle file name based on the workflow parameters
    network_basename = network_file.stem  # filename without extension
    pkl_pattern = f"active_no_parallel_p_0.02_max_step_{args.lambda_max_step}_hops_3_DE_DE_{dataset_stem}_{source_name}_to_{target_name}_adjacency_{network_basename}.pkl"
    plot_input_file = project_directory / pkl_pattern
    
    # Check if file exists, otherwise search for similar files
    if not plot_input_file.exists():
        print(f"[WARNING] Expected plot input file not found: {plot_input_file}")
        pkl_files = list(project_directory.glob("*.pkl"))
        if pkl_files:
            plot_input_file = pkl_files[-1]  # Use most recent pkl file
            print(f"[INFO] Using alternative pkl file: {plot_input_file}")
        else:
            raise FileNotFoundError(f"No .pkl files found in {project_directory}")
    
    make_precision_plot(
        fname=plot_input_file,
        max_rows=args.max_rows,
        source=source_name.capitalize(),
        target=target_name.capitalize(),
        out_dir=project_directory,
        out_name=args.plot_output,
    )
    
    # =========================================================================
    # 7. (OPTIONAL) EXPORT GRN FOR VISUALIZATION
    # =========================================================================
    # NOTE: NOT FUNCTIONAL YET (2025-09-29)
    # print("[INFO] Exporting network for visualization...")
    # regulon_path = project_directory / "allcells"
    # args = SimpleNamespace(
    #     regulons=str(regulon_path),
    #     species="human",
    #     name=None
    # )
    # create_network.create_net(args)
    
    print("[INFO] PreciCE workflow completed successfully!")


if __name__ == "__main__":
    main()
