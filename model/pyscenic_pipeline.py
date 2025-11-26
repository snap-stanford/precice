#!/usr/bin/env python
"""
PySCENIC end‑to‑end pipeline *for PreciCE*
================================================
A **pure‑Python** replacement for the original shell workflow that:

* optionally **sub‑selects cells** by a metadata field (default: `celltype`).
* runs GRN inference (GRNBOOST2) → motif pruning (pySCENIC) → regulon export.
* writes exactly the CSV edge‑list structure used by PreciCE’s weight‑learning
  step – e.g. `embryonic_stem_cells_meso_endo_Friedman_1_linearweights.csv` –
  including the left‑hand unnamed index column.

Usage (inside any notebook / script):

```python
from pyscenic_pipeline import run_pyscenic
run_pyscenic(
    adata      = my_adata,
    prefix     = "Friedman_1",  # becomes part of the file name ‑ see below
    species    = "human",
    metadata_key = "celltype",  # obs column to prompt on
    n_workers  = 16,
    seed       = 777,
    out_dir    = "./results"     # will be created if absent
)
```

The function will pause and show the distinct values in `adata.obs[metadata_key]`,
then prompt:
```
Enter comma‑separated cell types to include (or "all"): 
```
Press **Enter** for *all* cells. The selected subset is used for the GRN and
its label becomes part of the output filename::

    <celltypes>_<prefix>_linearweights.csv

So if you choose the three groups `embryonic_stem_cells,meso,endo` and keep the
`prefix='Friedman_1'`, the file written is::

    embryonic_stem_cells_meso_endo_Friedman_1_linearweights.csv
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path
from typing import Iterable, Sequence, Optional, List

import numpy as np
import pandas as pd

from scipy import sparse #Required for handling sparse input data (2025-05-29)

import learn_weights #For post-processing the GRN using regression (2025-06-26)
from types import SimpleNamespace #For creating path arguments isolated to one object (for learn_weights)

# ──────────────────────────────────────────────────────────────────────────────
#  PySCENIC / Arboreto imports (heavy)
# ──────────────────────────────────────────────────────────────────────────────
from arboreto.algo import grnboost2
from pyscenic.aucell import aucell
from pyscenic.cli.utils import load_signatures
from pyscenic.cli.utils import save_enriched_motifs
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase


from pyscenic.utils import load_motifs #Removing as test
from pyscenic.utils import modules_from_adjacencies

from dask.distributed import Client, LocalCluster #2025-05-30, for parallel processing of GRN-BOOST2
from threadpoolctl import threadpool_limits #force only the GRNBOOST2 phase to run with a single BLAS /OpenMP thread by wrapping the call in a threadpoolctl context manager. That keeps the rest of your Python session (e.g. SCVI training) free to use its original thread settings.

###New 2025-05-31
from pyscenic.prune import prune2df
###

# PySCENIC expects some global resources (TF list + motif database). These will be fetched from the precice/input/resources folder
DEFAULT_RESOURCE_DIR = Path.cwd().parent / "data" / "input" / "resources"

DEFAULT_TF_LIST = DEFAULT_RESOURCE_DIR / "TF_names_v_1.01_human.txt"
DEFAULT_MOTIF_DB = [
    DEFAULT_RESOURCE_DIR / "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
]
DEFAULT_MOTIF_ANNOTATIONS = (
    DEFAULT_RESOURCE_DIR / "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
)

# Safety check: warn if expected resource files are missing
for f in [DEFAULT_TF_LIST, *DEFAULT_MOTIF_DB, DEFAULT_MOTIF_ANNOTATIONS]:
    if not f.exists():
        print(f"[WARN] Resource file not found: {f}")
    else:
        print(f"[INFO] Found resource file: {f}")

# ──────────────────────────────────────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _prompt_cell_types(adata, metadata_key: str = "celltype") -> Optional[List[str]]:
    """Interactively ask the user to choose which groups in *adata.obs[metadata_key]* to keep.
    Returns **None** for all cells or a list of selected categories."""
    if metadata_key not in adata.obs:
        print(f"[WARN] metadata_key='{metadata_key}' not found in adata.obs – using all cells.")
        return None
    cats = sorted(adata.obs[metadata_key].astype(str).unique())
    print("\nAvailable categories in adata.obs['{0}'] (total {1}):".format(metadata_key, len(cats)))
    for c in cats:
        print(f"  • {c}")
    inp = input("\nEnter comma‑separated cell types to include (or 'all' for every cell): ").strip()
    if not inp or inp.lower() == "all":
        return None  # keep everything
    chosen = [x.strip() for x in inp.split(',') if x.strip() in cats]
    if not chosen:
        print("[WARN] No valid categories entered – using all cells.")
        return None
    print(f"→ Subsetting to {len(chosen)} categories: {', '.join(chosen)}")
    return chosen


def _subset_adata(adata, cell_types: Optional[Sequence[str]], metadata_key: str) -> "AnnData":
    if cell_types is None:
        return adata  # no subsetting
    mask = adata.obs[metadata_key].astype(str).isin(cell_types)
    return adata[mask].copy()

# ──────────────────────────────────────────────────────────────────────────────
#  Public API
# ──────────────────────────────────────────────────────────────────────────────
### NEW HELPER FUNCTION TO CONVERT OUTPUT REGULONS TO FLATTENED EDGE LIST TO ALLOW EXPORT TO EDGELIST (CORRECTS ERROR WHERE RAW ADJACENCIES WERE EXPORTED INSTEAD OF PRUNED NETWORK) 2025-06-25 ###
def regulons_to_edgelist(regulons):
    """
    Converts a regulons DataFrame from prune2df to an edge list DataFrame
    with columns: TF, target, importance.
    """
    edge_rows = []
    for (tf, motif_id), row in regulons.iterrows():
        target_tuples = row[('Enrichment', 'TargetGenes')]
        # Defensive: sometimes it's a list of tuples, sometimes list of strings
        if isinstance(target_tuples, list) and target_tuples and isinstance(target_tuples[0], tuple):
            for gene, weight in target_tuples:
                edge_rows.append({'TF': tf, 'target': gene, 'importance': weight})
        elif isinstance(target_tuples, list):
            for gene in target_tuples:
                edge_rows.append({'TF': tf, 'target': gene, 'importance': None})
        else:
            print(f"[WARN] Unexpected TargetGenes for TF={tf}: {type(target_tuples)}")
    return pd.DataFrame(edge_rows)

def run_pyscenic(
    *,
    adata,
    prefix: str,
    species: str = "human",
    metadata_key: str = "celltype",
    cell_types: Optional[Sequence[str]] = None,
    n_workers: int = 16,
    seed: int = 777,
    tf_list: Path | str = DEFAULT_TF_LIST,
    motif_dbs: Sequence[Path | str] = DEFAULT_MOTIF_DB,
    motif_annotations: Path | str = DEFAULT_MOTIF_ANNOTATIONS,
    out_dir: Path | str = "./results",
    interactive: bool = True,
):
    """Run the full GRN pipeline and write a PreciCE‑style edge list.

    Parameters
    ----------
    adata : AnnData | pd.DataFrame
        Raw counts matrix (*cells × genes*). If an AnnData, we assume genes in
        `adata.var_names` and cells in `adata.obs_names`.
    prefix : str
        Label that identifies this run (becomes part of the output file name).
    species : {"human", "mouse"}
        Used only for informative printout; *you* supply the correct TF list &
        motif databases for the species.
    metadata_key : str, default "celltype"
        Column in `adata.obs` that holds the user‑selectable categories.
    cell_types : list[str] | None
        Explicit list of categories to keep. If **None** and *interactive* is
        True, the user is prompted. If None and interactive is False, all cells
        are kept.
    interactive : bool, default True
        Whether to drop into an *input()* prompt when *cell_types* is None.
    out_dir : Path | str
        Directory in which *all* outputs are written.
    """

    from anndata import AnnData  # local import to avoid mandatory dependency if not needed

    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

# 1 ─────────────────────────  Subset
    if cell_types is None and interactive:
        cell_types = _prompt_cell_types(adata, metadata_key=metadata_key)
    adata_sub = _subset_adata(adata, cell_types, metadata_key)
    subset_tag = "_".join(cell_types) if cell_types else "allcells"

    print(f"[INFO] Expression matrix after subsetting: {adata_sub.shape[0]} cells × {adata_sub.shape[1]} genes")

    # ────────────────────────────────────────────────────────────────
    # Build expression matrix  (genes × cells)  and clean identifiers
    # ────────────────────────────────────────────────────────────────

    # ────────────────────────────────────────────────────────────────
    # Build expression matrix   (genes × cells)   and harmonise names
    # ────────────────────────────────────────────────────────────────
    if isinstance(adata_sub, AnnData):
        ex_matrix = adata_sub.to_df()          # cells × genes
    else:
        ex_matrix = adata_sub                  # assume DataFrame

    ex_matrix = ex_matrix.T                    # genes × cells

    # ① drop duplicate gene symbols, ② strip + upper-case
    ex_matrix = (
        ex_matrix
          .loc[~ex_matrix.index.duplicated(keep="first")]
          .rename(index=lambda s: s.strip().upper())
    )

    # Clean TF list the same way
    tf_names = (
        pd.read_csv(tf_list, header=None)[0]
          .str.strip()
          .str.upper()
          .drop_duplicates()
          .tolist()
    )

    #### INFO/FAIL-TEST BLOCK ###
    overlap = set(ex_matrix.index).intersection(tf_names)
    print("[INFO] first 5 genes :", list(ex_matrix.index[:5]))
    print("[INFO] first 5 TFs   :", tf_names[:5])
    print(f"[INFO] overlap size  : {len(overlap)}")
    if not overlap:
        raise RuntimeError("No overlap between gene index and TF list after cleaning.")
    print("[INFO] actual ex_matrix index dtype and values:", ex_matrix.index[:5])
    print("[INFO] gene_names[:5] before passing to grnboost2:", ex_matrix.index.tolist()[:5])
    print(set(map(type, tf_names)))
    print(repr(tf_names[:5]))
    print(repr(ex_matrix.index.tolist()[:5]))
    print("arboreto tf_names[:5]:", tf_names[:5])
    print("arboreto expression_matrix.index[:5]:", ex_matrix.index[:5])

    # 2 ─────────────────────────  GRN inference (GRNBOOST2)
    print(f"[INFO] Running GRNBOOST2 with {n_workers} workers …")

    with Client(LocalCluster(n_workers=n_workers, threads_per_worker=1)) as dask_client: #This function was written with ChatGPT 2025-05-29 and is intended to run GRN-BOOST on a set number of cores while capping number of threads to 1 per worker
        # cap BLAS / OpenMP threads to 1 **inside each worker**
        with threadpool_limits(limits=1):        
            adjacencies = grnboost2(
                ex_matrix.T,
                tf_names=tf_names,
                gene_names=ex_matrix.columns.tolist(),
                client_or_address=dask_client,
                seed=seed,
            )

    # 3 ─────────────────────────  Build and prune modules with motif info
    print("[INFO] Building modules from adjacencies …")
    
    modules = list(modules_from_adjacencies(
        adjacencies,
        ex_matrix.T,
        rho_mask_dropouts=True #Current pySCENIC default = False. Set to True to retain EARLIER pySCENIC default. True has always been used for PrecicE. True will only calculate correlations with genes that have non-zero expression values. May increase sensitivity for lowly expressed genes.
    ))

    dbs = [RankingDatabase(fname=str(db), name=db.stem) for db in motif_dbs]
    
    ### INFO BLOCK ###
    print(f"Modules before prune2df: {len(modules)}")
    ### END INFO BLOCK ###

### NOTE: THE FOLLOWING ARE DEFAULT SETTINGS ###
    regulons = prune2df( #replicates original values from the CLI pyscenic ctx command from first PreciCE release
        dbs, 
        modules, 
        str(motif_annotations),
        rank_threshold=5000,        # same defaults the CLI used
        auc_threshold=0.05,
        nes_threshold=3.0,
        motif_similarity_fdr=0.001,
        filter_for_annotation=True,
        orthologuous_identity_threshold=0.0,
        module_chunksize=100,
        weighted_recovery=False
    )
### END DEFAULT SETTING BLOCK ###

### NOTE: THE FOLLOWING ARE PERMISSIVE SETTINGS FOR REGULON PRUNING -> MORE/LARGER REGULONS ###
##Permissive settings
    # regulons = prune2df( 
    #     dbs, 
    #     modules, 
    #     str(motif_annotations),
    #     rank_threshold=6000,        
    #     auc_threshold=0.03,
    #     nes_threshold=2.5,
    #     motif_similarity_fdr=0.01,
    #     filter_for_annotation=True,
    #     orthologuous_identity_threshold=0.0,
    #     module_chunksize=100,
    #     weighted_recovery=False
    # )
### END PERMISSIVE SETTING BLOCK ###


    ### INFO BLOCK ###
    print("regulons.columns:", regulons.columns)
    print("regulons.head():", regulons.head())
    print("regulons.shape:", regulons.shape)
    print("Regulons shape:", regulons.shape)
    ### END INFO BLOCK ###

    # Now create the edge list from pruned regulons
    edge_df = regulons_to_edgelist(regulons)


    # 4 ─────────────────────────  Write outputs expected by PreciCE
    adj_path = out_dir / f"{prefix}_{subset_tag}_adjacencies.csv"
    adjacencies.to_csv(adj_path, index=False)
    print(f"[DONE] Raw adjacencies → {adj_path}")

    # Edge list in PreciCE *linearweights* format (with unnamed index column) (CORRECTED 2025-06-25)
    lw_path = out_dir / f"{prefix}_{subset_tag}_linearweights.csv"
    edge_df.to_csv(lw_path, index=False)
    print(f"[DONE] PreciCE edge list (from regulons) → {lw_path}")

    # Optionally export regulons in GMT + per‑TF text files
    reg_dir = out_dir / subset_tag
    reg_dir.mkdir(exist_ok=True)
    ### INFO BLOCK: Inspect regulons before for-loop ###
    import sys
    import inspect

    print("\n" + "="*40)
    print("INFO: Before writing regulon files")
    print(f"Python version: {sys.version}")
    try:
        import pyscenic
        print("pySCENIC version:", pyscenic.__version__)
        import ctxcore
        print("ctxcore version:", ctxcore.__version__)
    except Exception as e:
        print("pyscenic/ctxcore version check failed:", e)

    print("Type of regulons:", type(regulons))
    if hasattr(regulons, 'columns'):
        print("regulons columns:", regulons.columns)
        print("First row (iloc[0]):", regulons.iloc[0])
        print("type of first row:", type(regulons.iloc[0]))
        if 'regulon' in regulons.columns:
            print("type of regulons['regulon'][0]:", type(regulons['regulon'].iloc[0]))
    elif isinstance(regulons, list):
        print("Length of regulons:", len(regulons))
        if len(regulons) > 0:
            print("Type of regulons[0]:", type(regulons[0]))
            print("regulons[0]:", regulons[0])
            # If it's a tuple, print tuple element types
            if isinstance(regulons[0], tuple):
                print("Types in regulons[0] tuple:", [type(x) for x in regulons[0]])
    else:
        print("Unknown regulons type.")

    # Print head if DataFrame
    if hasattr(regulons, 'head'):
        print("First few rows of regulons:\n", regulons.head())

    # Print contents if list
    if isinstance(regulons, list):
        print("First few regulons (list):\n", regulons[:3])

    print("="*40 + "\n")
    ### END INFO BLOCK ###

### INFO BLOCK for "for reg in regulons" replacement ###
    print("\n[INFO] About to write regulons from DataFrame")
    print("regulons index (first 5):", list(regulons.index[:5]))
    print("regulons columns:", regulons.columns)
    print("Number of regulons:", len(regulons))

    bad_entries = 0

    for i, ((tf, motif_id), row) in enumerate(regulons.iterrows()):
        target_tuples = row[('Enrichment', 'TargetGenes')]
        # Info: show what we're seeing in the first row or two
        if i < 2:
            print(f"[INFO] Row {i}: TF={tf}, motif_id={motif_id}")
            print("  TargetGenes type:", type(target_tuples))
            if isinstance(target_tuples, list) and target_tuples:
                print("  First 3 targets:", target_tuples[:3])
            else:
                print("  TargetGenes value:", target_tuples)
        # Defensive extraction
        if isinstance(target_tuples, list) and target_tuples and isinstance(target_tuples[0], tuple):
            targets = sorted(gene for gene, weight in target_tuples)
        elif isinstance(target_tuples, list):
            targets = [str(g) for g in target_tuples]
            print(f"[WARN] Row {i}: TargetGenes list does not contain tuples for TF={tf}. Using as strings.")
            bad_entries += 1
        else:
            print(f"[ERROR] Row {i}: Unexpected TargetGenes type ({type(target_tuples)}), TF={tf}")
            targets = []
            bad_entries += 1
        fname = f"{tf}__{motif_id}.txt"
        try:
            (reg_dir / fname).write_text("\n".join(targets))
        except Exception as e:
            print(f"[ERROR] Failed to write file for TF={tf}, motif_id={motif_id}: {e}")

    print(f"[DONE] Regulon directory → {reg_dir}")
    if bad_entries:
        print(f"[WARN] {bad_entries} entries had unexpected TargetGenes values. Check previous warnings above.")
### END INFO BLOCK ###

    gmt_path = out_dir / f"{prefix}_{subset_tag}_regulons.gmt"
    save_enriched_motifs(regulons, str(gmt_path))
    print(f"[DONE] Regulons GMT → {gmt_path}")

    return lw_path  # path to the main CSV

__all__ = ["run_pyscenic"]
