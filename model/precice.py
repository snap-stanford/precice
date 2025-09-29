import os
import subprocess
from pathlib import Path

import scanpy as sc
import scvi
import pandas as pd
import numpy as np

from utils import mito_qc, keep_variable_genes, remove_mito_ribo_genes
from pyscenic_pipeline import run_pyscenic as _run_pyscenic

"""
This class orchestrates the PreciCE workflow, including preprocessing,
SCVI batch correction, differential expression, PySCENIC GRN inference,
and downstream network construction and perturbation ranking.
"""

class precice:
    def __init__(
        self,
        dir:str = './',
        path: str | None = None,
        adata=None,
        batch_correct: bool = False,
        cell_filter: bool = False,
        gene_filter: bool = True,
        all_network_paths_file: str = '../data/input/networks/all_network_paths.csv',
        species: str = 'human',
        mito_percent_max: float = 5,   # ← Added 2025-08-21 to enable passing max allowed mitochondrial reads in precice()
    ):
        # Basic setup
        self.name = Path(path).stem if path else 'run'
        self.path = path
        self.dir = dir
        self.adata = None
        self.tf_list = None
        self.f_motif_path = None
        self.f_db_name = None
        self.cell_filter = cell_filter
        self.gene_filter = gene_filter
        self.batch_correct = batch_correct
        self.mito_percent_max = mito_percent_max   # ← Added 2025-08-21 to enable passing max allowed mitochondrial reads in precice()
        self.source = None
        self.target = None
        self.corrected_subset = None
        self.trans_subset = None
        self.DE = {}
        self.DE_filenames = {}

        # Load network path mappings
        net_paths = pd.read_csv(all_network_paths_file, index_col=0)
        net_paths = net_paths[net_paths['species'] == species]
        self.all_network_paths = net_paths.set_index('cell_type')['path'].to_dict()

        # Load or preprocess AnnData
        if adata is not None:
            self.adata = adata
        if self.batch_correct:
            self.processed_name = os.path.join(self.dir, self.name + '_scvi_processed.h5ad')
        else:
            self.processed_name = os.path.join(self.dir, self.name + '_processed.h5ad')
        if os.path.isfile(self.processed_name):
            print('Loading preprocessed data.')
            self.adata = sc.read_h5ad(self.processed_name)
        else:
            if self.adata is None and self.path:
                print('Loading raw data.')
                self.adata = sc.read_h5ad(self.path)
            self.preprocess()
            #self.adata.write_h5ad(self.processed_name) #No longer necessary to output this processed file because I already output two specific ones below

    def preprocess(self):
        """
        Standard Scanpy preprocessing: filtering, normalization, log1p,
        high-variable gene selection, and mito/ribo gene removal.
        """
        print(f"[INFO] Start: {self.adata.n_obs} cells × {self.adata.n_vars} genes")
        if self.cell_filter:
            sc.pp.filter_cells(self.adata, min_genes=200)
            print(f"[INFO] After filter_cells(min_genes=200): {self.adata.n_obs} cells")
        if self.gene_filter:
            sc.pp.filter_genes(self.adata, min_cells=3)
            print(f"[INFO] After filter_genes(min_cells=3): {self.adata.n_vars} genes")
        if self.cell_filter:
            n_before = self.adata.n_obs
            self.adata = mito_qc(self.adata, mito_percent_max=self.mito_percent_max) #Edited 2025-08-21 for mitochondrial read specification
            print(f"[INFO] After mito_qc: {self.adata.n_obs} cells (removed {n_before - self.adata.n_obs})")
        self.adata.layers['counts'] = self.adata.X.copy()
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        print(f"[INFO] After normalization+log1p: {self.adata.n_obs} cells × {self.adata.n_vars} genes")
        if self.gene_filter:
            self.adata = remove_mito_ribo_genes(self.adata)
            print(f"[INFO] After remove_mito_ribo_genes: {self.adata.n_obs} cells × {self.adata.n_vars} genes")
            
            # Save the PySCENIC-ready file (no highly-variable-gene filter)
            scenic_fp = os.path.join(
                self.dir, f"{self.name}_processed_for_pyscenic.h5ad"
            )
            self.adata.write_h5ad(scenic_fp)

            # Now apply highly-variable-gene filtering for DE / scVI            
            self.adata = keep_variable_genes(self.adata) 
            print(f"[INFO] After keep_variable_genes: {self.adata.n_obs} cells × {self.adata.n_vars} genes")
            # Write the DE-ready file
            de_fp = os.path.join(
                self.dir, f"{self.name}_processed_for_DE.h5ad"
            )
            self.adata.write_h5ad(de_fp)


    def set_up_scvi(self, batch_key: str):
        """
        Initialize and train an scVI model for batch correction.
        """
        #scvi.model.SCVI.setup_anndata(self.adata, layer='counts', batch_key=batch_key) #Removed after dependency updates to be compatible with other version of scvi
        self.adata = self.adata.copy()
        scvi.model.SCVI.setup_anndata(self.adata, layer='counts', batch_key=batch_key)
        self.scvi_model = scvi.model.SCVI(self.adata, gene_likelihood='nb')
        self.scvi_model.train(
            check_val_every_n_epoch=1,
            max_epochs=400,
            early_stopping=True,
            early_stopping_patience=20,
            early_stopping_monitor='elbo_validation',
        )

    def scvi_plot_setup(self):
        """
        Embed corrected latent space for visualization.
        """
        latent = self.scvi_model.get_latent_representation()
        self.adata.obsm['X_scVI'] = latent
        sc.pp.neighbors(self.adata, use_rep='X_scVI')
        sc.tl.umap(self.adata)

    def batch_effect_correction(self, batch_key: str='batch', batches=None):
        """
        Perform batch correction using scGen for specified batches.
        """
        if batches is None:
            batches = list(self.adata.obs[batch_key].unique())
        subset = self.adata[self.adata.obs[batch_key].isin(batches)].copy()
        import scgen
        scgen.SCGEN.setup_anndata(subset, batch_key=batch_key)
        model = scgen.SCGEN(subset)
        model.train(
            max_epochs=100, batch_size=32,
            early_stopping=True, early_stopping_patience=25
        )
        self.corrected_subset = model.batch_removal()
        return self.corrected_subset

    def get_DE(
        self,
        source_idx: np.ndarray | None = None,
        target_idx: np.ndarray | None = None,
        target_name: str = "target", 
        source_name: str = "source", 
        precomputed_DE: str | None = None,
        cell_type_label: str = "label", 
        batch_key: str | None = "batch", 
        quick_train_epochs: int = 50,
        *,                         # -------- keyword-only after this
        de_method: str = "scvi",   # "scvi"  or "rank_genes_groups"
        rgg_method: str = "wilcoxon",
    ):
        """
        Differential expression between *source* and *target* groups.

        Parameters
        ----------
        de_method : {"scvi", "rank_genes_groups"}, default "scvi"
            • "scvi" – use a (possibly quick-trained) scVI generative model.\n
            • "rank_genes_groups" – use Scanpy’s rank_genes_groups (no SCVI).
        rgg_method : str, default "wilcoxon"
            Statistical test passed to ``scanpy.tl.rank_genes_groups`` when
            *de_method="rank_genes_groups"*.
        quick_train_epochs : int, default 50
            If no SCVI model exists and *de_method="scvi"*, train a light
            model for this many epochs.
        """
        trans_name = f"{source_name}_to_{target_name}"

        # ------------------------------------------------------------------
        # 0)  Return cached / pre-computed table
        # ------------------------------------------------------------------
        if precomputed_DE:
            self.DE[trans_name] = pd.read_csv(precomputed_DE, index_col=0)
            return self.DE[trans_name]

        # ------------------------------------------------------------------
        # 1)  Build source / target index masks
        # ------------------------------------------------------------------
        if source_idx is None:
            source_idx = np.where(
                self.adata.obs[cell_type_label] == source_name
            )[0]
        if target_idx is None:
            target_idx = np.where(
                self.adata.obs[cell_type_label] == target_name
            )[0]

        # ------------------------------------------------------------------
        # 2)  SCVI generative-model DE
        # ------------------------------------------------------------------
        if de_method == "scvi":
            if getattr(self, "scvi_model", None) is None:
                print(
                    f"[INFO] de_method='scvi' but no SCVI model found – "
                    f"training a quick model ({quick_train_epochs} epochs)."
                )
                self.adata = self.adata.copy()  # avoid AnnData view error
                scvi.model.SCVI.setup_anndata(
                    self.adata, layer="counts", batch_key=batch_key
                )
                self.scvi_model = scvi.model.SCVI(
                    self.adata, gene_likelihood="nb"
                )
                self.scvi_model.train(
                    max_epochs=quick_train_epochs,
                    early_stopping=True,
                    early_stopping_patience=10,
                    check_val_every_n_epoch=1,
                    early_stopping_monitor="elbo_validation",
                )

            de_raw = self.scvi_model.differential_expression(
                idx1=target_idx, idx2=source_idx
            )
            self.DE[trans_name] = self.convert_de_to_seurat_format(de_raw)
            return self.DE[trans_name]

        # ------------------------------------------------------------------
        # 3)  Scanpy / Seurat-style DE
        # ------------------------------------------------------------------
        elif de_method == "rank_genes_groups":
            import scanpy as sc

            # temporary subset with just the two groups
            a_tmp = self.adata[[*source_idx, *target_idx], :].copy()
            a_tmp.obs["tmp_group"] = [
                source_name if i < len(source_idx) else target_name
                for i in range(a_tmp.n_obs)
            ]

            sc.tl.rank_genes_groups(
                a_tmp,
                groupby="tmp_group",
                groups=[target_name],
                reference=source_name,
                method=rgg_method,
            )
            df = sc.get.rank_genes_groups_df(a_tmp, group=target_name)
            df = df.rename(
                columns={
                    "names": "gene_name",
                    "logfoldchanges": "avg_log2FC",
                    "pvals": "p_val",
                    "pvals_adj": "p_val_adj",
                }
            )
            df["pct.1"] = -1
            df["pct.2"] = -1
            self.DE[trans_name] = df
            return df

        # ------------------------------------------------------------------
        else:
            raise ValueError(
                "de_method must be 'scvi' or 'rank_genes_groups'; "
                f"got '{de_method}'."
            )

    def write_DE_files(self, DE_dir: str='./'):
        """
        Write DE results to CSV files.
        """
        for name, df in self.DE.items():
            out = os.path.join(DE_dir, f'DE_{self.name}_{name}.csv')
            df.to_csv(out)
            self.DE_filenames[name] = out

    def convert_de_to_seurat_format(self, de_df):
        """
        Specifically designed for scVI output
        """

        de_df = de_df.reset_index()
        de_df['p_val_adj'] = 1.0
        de_df = de_df[de_df['is_de_fdr_0.05']==True]
        de_df = de_df[np.logical_or(de_df['non_zeros_proportion1']>0.5, de_df['non_zeros_proportion2']>0.5)]
        de_df = de_df[np.logical_or(de_df['raw_mean1']>1, de_df['raw_mean2']>1)]

        ## Since significance test was performed above just set value to 0.01
        # so it isn't filtered by PreciCE
        de_df['p_val_adj'] = 0.01
        de_df['p_val'] = 0.01
        de_df['pct.1'] = -1.0
        de_df['pct.2'] = -1.0

        de_df = de_df.rename(columns={'index':'gene_name',
                              'lfc_mean': 'avg_log2FC'})
        de_df = de_df.loc[:, ['gene_name', 'p_val',
                              'avg_log2FC', 'pct.1','pct.2',
                              'p_val_adj']]
        return de_df

    def run_pyscenic(
        self,
        species: str = 'human',
        metadata_key: str = 'label',
        n_workers: int = 16,
        seed: int = 777,
        out_dir: str | Path = '/home/jensm/Programs/miniconda3/envs/pyscenic_env/precice/data/output'
    ):
        """
        Infer a gene regulatory network via PySCENIC in pure-Python.
        """
        _run_pyscenic(
            adata=self.adata,
            prefix=self.name,
            species=species,
            metadata_key=metadata_key,
            n_workers=n_workers,
            seed=seed,
            out_dir=out_dir,
            interactive=True,
        )

    def get_network(
        self,
        cell_type: str | None = None,
        network_path: str | Path | None = None,
    ):
        """
        Register the GRN file to be used downstream.

        Parameters
        ----------
        cell_type
            Key into ``self.all_network_paths`` (old behaviour).
            Ignored if *network_path* is supplied.
        network_path
            Full path to a .csv network file.
            If provided, we skip the lookup table.
        """
        if network_path is not None:
            # user gave us an explicit file → normalise & store it
            self.network_path = os.fspath(network_path)
        elif cell_type is not None:
            # fall back to the mapping that was already in place
            path = self.all_network_paths[cell_type]
            # honour the project folder if self.dir exists
            base = getattr(self, "dir", "../data/input/networks")
            self.network_path = os.path.join(base, path)
        else:
            raise ValueError("Provide either 'network_path' or 'cell_type'.")

        print("Network path loaded:", self.network_path)

    def run_precice(
        self,
        species: str,
        network_path: str,
        DE_path: str,
        python_path: str,
        *,
        out_dir: str | Path | None = None,
        lambda_max_step: int  | None = None,    # Allows passing lambda_max_step in run_precice() (2025-07-09)
        lambda_step_size: float | None = None,  # Allows passing lambda_step_size in run_precice() (2025-07-09)
        pos_only: bool = False,
        remove_TFs: list[str] | None = None,
    ):
        """
        Run the core PRECICE solver with given network and DE data.
        """
        if out_dir is None:
            out_dir = getattr(self, "dir", ".")      # default to project_directory if present

        cmd = [python_path,
               './run_precice.py',
               '--p_thresh','2e-2','--hops','3',
               '--species',species,
               '--adjacency',network_path,
               '--DE_data',DE_path,
               '--out_dir', os.fspath(out_dir)]      # ← NEW
        ### ---- forward the optional lambda knobs: -------------
        if lambda_max_step  is not None:
            cmd.extend(['--max_step',  str(lambda_max_step)])
        if lambda_step_size is not None:
            cmd.extend(['--step_size', str(lambda_step_size)])
        ###
        if pos_only: cmd.append('--pos_pert')
        if remove_TFs: cmd.extend(['--remove_TFs', *remove_TFs])
        print('Running PRECICE:', ' '.join(cmd))
        subprocess.run(cmd)
