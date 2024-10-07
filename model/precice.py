import scanpy as sc
import scvi
import pandas as pd
from utils import remove_noise, mito_qc, keep_variable_genes, \
    remove_mito_ribo_genes
import numpy as np
#import scgen
from pyscenic_utils import create_loom, get_non_unique_IDs, save_loom,\
create_target_genes, save_regulons
import subprocess
import os

"""
This class is used to preprocess the data and run PRECICE. It is also used 
to store the data and the results of the analysis. 

"""

class precice():

    def __init__(self, dir='./', path=None, adata=None, batch_correct=False, 
                 cell_filter=False, gene_filter=True, 
                 all_network_paths_file='../data/all_network_paths.csv',
                 species='human'):

        self.name = path.split('/')[-1].split('.')[0]
        self.path = path
        self.dir = dir
        self.adata = None
        self.tf_list = None
        self.f_motif_path = None
        self.f_db_name = None
        
        self.cell_filter = cell_filter
        self.gene_filter = gene_filter
        self.batch_correct = batch_correct

        self.source = None
        self.target = None
        self.corrected_subset = None
        self.trans_subset = None
        self.DE = {}
        self.DE_filenames = {}
        
        # Set up network paths
        net_paths = pd.read_csv(all_network_paths_file, index_col=0)
        net_paths = net_paths[net_paths['species'] == species]
        self.all_network_paths = net_paths.set_index('cell_type').loc[:, 'path'].to_dict()

        # check if processed file exists
        if adata is not None:
            self.adata = adata
            if self.path is None:
                raise ValueError('Please provide path to raw data file')

        if batch_correct:
            self.scvi_model = None
            self.processed_name = os.path.join(self.dir,self.name + '_scvi_processed.h5ad')
        else:
            self.processed_name = os.path.join(self.dir,self.name + '_processed.h5ad')
        if os.path.isfile(self.processed_name): 
            print('Loading preprocessed data...')
            self.adata = sc.read_h5ad(self.processed_name)
        else:
            print('Loading raw data...')
            self.adata = sc.read_h5ad(self.path)
            self.preprocess()
            self.adata.write_h5ad(self.processed_name)
            
    def preprocess(self):
        """
        Preprocess the data using standard scanpy functions
        """
        ## TODO Check if var.index is gene names and is unique
        
        if self.cell_filter:
            sc.pp.filter_cells(self.adata, min_genes=200)
        
        if self.gene_filter:
            sc.pp.filter_genes(self.adata, min_cells=3)
            
        if self.cell_filter:
            self.adata = mito_qc(self.adata)

        #scvi.data.poisson_gene_selection(self.adata, n_top_genes=n_top_genes)
        #self.adata = self.adata[:, self.adata.var["highly_variable"]] 
        #self.adata.layers['counts'] = self.adata.X.copy()
            
        
        self.adata.layers['counts'] = self.adata.X.copy()
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
            
        if self.gene_filter:
            self.adata = keep_variable_genes(self.adata)
            self.adata = remove_mito_ribo_genes(self.adata)
            

    def set_up_scvi(self, batch_key):
        """
        Set up scvi model if performing batch correction using scVI
        """
        self.adata = self.adata.copy()
        scvi.model.SCVI.setup_anndata(
            self.adata, layer="counts", batch_key=batch_key
        ) 
        
        self.scvi_model = scvi.model.SCVI(
            self.adata, gene_likelihood="nb"
        ) 
        
        self.scvi_model.train(
            check_val_every_n_epoch=1,
            max_epochs=400,
            early_stopping=True,
            early_stopping_patience=20,
            early_stopping_monitor="elbo_validation",
        )
        
    def scvi_plot_setup(self):
        
        SCVI_LATENT_KEY = "X_scVI"
        latent = self.scvi_model.get_latent_representation()
        self.adata.obsm[SCVI_LATENT_KEY] = latent
        
        sc.pp.neighbors(self.adata, use_rep=SCVI_LATENT_KEY)
        sc.tl.umap(self.adata)

    def batch_effect_correction(self,  batch_key='batch', batches=None):
        """
        Batch effect correction using scGen
        """

        if batches is None:
            batches = self.adata.obs[batch_key].unique()
        subset = self.adata[self.adata.obs[batch_key].isin(batches)]
        subset = subset.copy()

        scgen.SCGEN.setup_anndata(subset, batch_key=batch_key)
        model = scgen.SCGEN(subset)

        model.train(
            max_epochs=100,
            batch_size=32,
            early_stopping=True,
            early_stopping_patience=25,
        )

        model.adata.obs = model.adata.obs.reset_index()
        self.corrected_subset = model.batch_removal()
        return self.corrected_subset

    def set_up_pyscenic(self, species, adata=None):
        """
        Set up PySCENIC for inferring gene regulatory networks
        """

        if adata is None:
            adata = self.adata
        reformat_data = adata.to_df().T.reset_index().rename(
            columns={'index': 'gene_name'})
        
        print("Writing .txt file to disk for PySCENIC")
        reformat_data.to_csv(self.name +'.txt', sep='\t', index=False)

        if species == 'mouse':
            # Load the transcription factor list
            self.tf_file = '../Data/TF_names/mouse_tf_gene_names_lowercase.txt'

            # motif database:
            self.f_motif_path = \
                '../Data/motifs/motifs-v9-nr.mgi-m0.001-o0.0.tbl'

            # ranking database
            self.f_db_name = '../Data/motifs/v2_database/mm10__refseq' \
                         '-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'

        elif species == 'human':
            self.tf_file = '../Data/TF_names/TF_names_v_1.01_human.txt'
            self.f_motif_path = '../Data/motifs/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
            self.f_db_name = '../Data/motifs/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'

        tfdb = pd.read_csv(self.tf_file, sep='\t', header=None)
        self.tf_list = tfdb.loc[:, 0].values

        self.loom_paths = {}
        self.loom_paths[self.name] = create_loom(self.name+'.txt')

    def set_transition(self, source_idx=None, target_idx=None,
                       label_map=None, index=True, colname='celltype',
                       source_name=None, target_name=None):
        """
        Set the source and target cells for the transition
        """

        label_map = {x:'Source' for x in source_idx}
        label_map.update({x:'Target' for x in target_idx})
        label_map_inv = {v: k for k, v in label_map.items()}

        if self.corrected_subset is not None:
            adata = self.corrected_subset
        else:
            adata = self.adata

        if index==True:
            adata.obs['Type'] = adata.obs.index.map(label_map)
        else:
            adata.obs['Type'] = adata.obs[colname].map(label_map)
        self.trans_subset = adata[~adata.obs['Type'].isna()]

        if source_name is None:
            self.source_name = label_map_inv['Source'].replace(' ', '_')
        else:
            self.source_name = source_name
            
        if target_name is None:
            self.target_name = label_map_inv['Target'].replace(' ', '_')
        else:
            self.target_name =target_name
        self.transition = self.name + '_' + self.source_name + '_' + self.target_name

    def save_seurat(self):
        """
        Save the data in a format that can be easily loaded into Seurat for
        differential expression analysis
        """

        self.trans_subset.write_h5ad(self.name + '_' + self.source_name + '_' +
                                self.target_name + '.h5ad')

        try:
            self.trans_subset.X = self.trans_subset.X.toarray()
        except:
            pass

        pd.DataFrame(data=self.trans_subset.X,
                        index=self.trans_subset.obs_names,
                        columns=self.trans_subset.var_names).to_csv(
            os.path.join(self.dir, self.transition + '.csv'), sep='\t')

        ## Save the data for Seurat
        self.trans_subset.write_csvs(os.path.join(self.dir, self.transition))

    def get_DE(self, source_idx=None, target_idx=None, 
               target_name='target', source_name='source', 
               precomputed_DE=None):
        """
        Compute differential expression between source and target cells using scVI
        """

        trans_name = source_name + '_to_' + target_name
        if precomputed_DE is not None:
            self.DE[trans_name] = pd.read_csv(DE, index_col=0)

        else:
            # Compute differential expression using scVI
            self.DE[trans_name] = self.scvi_model.differential_expression(idx1=target_idx, 
                                                                          idx2=source_idx)
            self.DE[trans_name] = self.convert_de_to_seurat_format(self.DE[trans_name])
            
    def write_DE_files(self, DE_dir='./'):
        """
        Write the differential expression results from scVI to a file
        """
        
        for DE_name, DE_file in self.DE.items():
            out_name = DE_dir + 'DE_' + self.name + '_' + DE_name +'.csv'
            print(out_name)
            DE_file.to_csv(out_name)
            self.DE_filenames[DE_name] = out_name
        
    
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

    def run_pyscenic(self, arboreto_path=None, pyscenic_path=None):
        """
        Run PySCENIC to infer gene regulatory networks
        """

        ## Step 1
        # Run GRNBOOST2
        ## Check if adjacency matrix exists
        ## if not then run arboreto and compute regulons
        
        print("Run these commands in the command line")

        for name, loom_path in self.loom_paths.items():

            if arboreto_path is None:
                arboreto_path = 'python3.7 ../../bin/arboreto_with_multiprocessing.py'

            if pyscenic_path is None:
                pyscenic_path = 'pyscenic'

            adjacency_file = 'adjacencies_' + name + '.csv'

            print(f"!{arboreto_path}\
                  {loom_path}\
                  {self.tf_file}\
                  --method grnboost2\
                  -o {adjacency_file}\
                  --num_workers 8\
                  --seed 777 \n")

        ## Step 2
        # Identify regulons
        for name, loom_path in self.loom_paths.items():
            adjacency_file = 'adjacencies_' + name + '.csv'
            regulons_path = 'regulons_' + name + '.csv'

            print(f"\n !{pyscenic_path} ctx ./{adjacency_file} \
            {self.f_db_name} \
            --annotations_fname {self.f_motif_path} \
            --expression_mtx_fname {loom_path} \
            --output ./{regulons_path} \
            --mask_dropouts \
            --num_workers 100")


    def save_pyscenic_output(self):

        # Save regulons output from PySCENIC analysis in appropriate format
        save_regulons(self.name)
        
    def get_network(self, cell_type):
        """
        Get the network path for precomputed network for a specific cell type
        """
        self.network_path = os.path.join('./data/networks',self.all_network_paths[cell_type])
        print('Network path loaded')


    def create_network(self):
        ## NOT implemented yet for website
        ## Use output from PySCENIC to create network
        # Sample usage: python create_network.py
        # --regulon ../Data/regulons_all/pijuan_sala_meso
        # --species mouse
        # --name pijuan_sala_meso
        pass

    def learn_weights(self):
        ## NOT implemented yet for website
        ## Use the expression matrix to learn the edge weights
        ## python learn_weights_updated.py
        # --data_path ../Notebooks/pijuan_sala_HVG.h5ad
        # --graph_name ../Data/transcription_networks/G_all_edges_pijuan_sala_meso
        # --out_dir ./
        pass

    def run_precice(self, species,
                      network_path, DE_path, python_path,
                      pos_only=False, remove_TFs=None):
        """
        Run PRECICE to infer ranking of genes to perturb to achieve desired
        transition
        """
        
        command_list = [python_path,
                        "../model/run_precice.py",
                                "--p_thresh", "2e-2",\
                                "--hops", "3",\
                                "--species", species,\
                                "--adjacency", network_path,\
                                "--DE_data", DE_path]
        
        if pos_only:
            command_list.append("--pos_pert")
        if remove_TFs is not None:
            to_add = ["--remove_TFs"] + remove_TFs
            command_list.extend(to_add)

            
        print("Running command:", ' '.join(command_list))
        result = subprocess.run(command_list)

