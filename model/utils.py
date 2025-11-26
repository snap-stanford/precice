import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

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


def mito_qc(adata, mito_percent_max=5):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate mito genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

    adata = adata[adata.obs.pct_counts_mt < mito_percent_max, :]
    return adata[adata.obs.pct_counts_mt < mito_percent_max, :].copy() #Modified 2025-08-21

def keep_variable_genes(adata):
    min_mean = 0.0125    #Default 0.0125 (used by Yusuf). Lower to include more lowly expressed genes (e.g., TFs). Try lowering to 0.005 to retain more TFs. I now lowered to 0.003 to make very permissive. but not as low as 0.001, which could include too much noise. Default: 0.0125
    max_mean = 3
    min_disp = 0.5      #Default 0.5 (used by Yusuf). 0.3 was my "Permissive". Now setting to very permissive (0.2)
    
    sc.pp.highly_variable_genes(adata, 
        min_mean=min_mean, 
        max_mean=max_mean,
        min_disp=min_disp 
        )
    adata = adata[:, adata.var.highly_variable].copy() #Added 2025-08-21
    return adata #Added 2025-08-21
    
    print("min_mean:", min_mean, "min_disp:", min_disp)
    print("Number of genes with mean >= min_mean:", np.sum(adata.var['means'] >= min_mean))
    print("Number of genes with dispersion >= min_disp:", np.sum(adata.var['dispersions'] >= min_disp))
    print("INFO: number of retained genes: ", sum(adata.var.highly_variable)) #Debug to show how many genes were retained after filtering

    sc.pl.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]
    return adata

def remove_mito_ribo_genes(adata):
    mito_genes = adata.var_names.str.startswith('MT-')
    ribo1_genes = adata.var_names.str.startswith('RPS')
    ribo2_genes = adata.var_names.str.startswith('RPL')

    remove = np.add(mito_genes, ribo1_genes)
    remove = np.add(remove, ribo2_genes)

    keep = np.invert(remove)

    return adata[:, keep].copy() #Modified 2025-08-21
