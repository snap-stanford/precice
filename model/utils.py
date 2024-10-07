import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import sys
import pandas as pd
import scanpy as sc
import loompy as lp
import numpy as np
import pickle
#import pyscenic
import networkx as nx
import os

def invert_dict(map_):
    return {v: k for k, v in map_.items()}

def get_gene_names(map_, list_):
    gene_names = []
    for l in list_:
        try:
            gene_names.append(invert_dict(map_)[l])
        except:
            gene_names.append(None)
    return gene_names

def ent_wcc_plot(ents, wccs, title):
    plt.title(title)
    #plt.plot(ents)
    plt.plot(wccs)
    plt.xlabel('Number of nodes removed')
    plt.ylabel('Number of components')
    plt.show()

def float_to_hsv(val, norm, cmap):
    m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mm = m.to_rgba(val)
    return colorsys.rgb_to_hsv(mm[0], mm[1], mm[2])


def remove_noise(adata, cell_filter=True, gene_filter=True):
    if cell_filter:
        sc.pp.filter_cells(adata, min_genes=200)
    if gene_filter:
        sc.pp.filter_genes(adata, min_cells=3)


def mito_qc(adata):
    adata.var['mt'] = adata.var_names.str.startswith(
        'MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

    # adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    return adata


def keep_variable_genes(adata):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                                min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]
    return adata
    # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # sc.pp.scale(adata, max_value=10)


def remove_mito_ribo_genes(adata):
    mito_genes = adata.var_names.str.startswith('MT-')
    ribo1_genes = adata.var_names.str.startswith('RPS')
    ribo2_genes = adata.var_names.str.startswith('RPL')

    remove = np.add(mito_genes, ribo1_genes)
    remove = np.add(remove, ribo2_genes)

    keep = np.invert(remove)

    adata = adata[:, keep]
    return adata
