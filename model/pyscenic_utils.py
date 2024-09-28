import pandas as pd
import scanpy as sc
import loompy as lp
import numpy as np
import os


# Helper functions
# Create loom file for genie
def create_loom(filename):
    adata = sc.read_text(filename, delimiter = '\t', first_column_names = True)
    loom_path = filename.split('.txt')[0 ] +'.loom'
    row_attrs_stem = {"Gene": np.array(adata.obs.index)}
    col_attrs_stem = {"CellID": np.array(adata.var.index)}
    adata.obs_names_make_unique()
    lp.create(loom_path, adata.X, row_attrs_stem, col_attrs_stem)
    return loom_path


def get_non_unique_IDs():
    return np.load('non_unique_IDs.txt.npy')


# Save the output as loom files
def save_loom(loom_file, name):
    reg_mtx_file = './reg_mtx_' + name
    auc_mtx_file = './auc_mtx_' + name

    lf = lp.connect(loom_file, mode='r+', validate=False)
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    reg_mtx = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)
    lf.close()

    auc_mtx.to_csv(auc_mtx_file, sep='\t')
    reg_mtx.to_csv(reg_mtx_file, sep='\t')
    
    return auc_mtx, reg_mtx


def create_target_genes(tflist):
    fixed = []
    for tf in tflist:
        genelist = tf[2:-2].split('), (')
        temp = []
        for gene in genelist:
            symb = gene.split(',')[0][1:-1]
            temp.append(symb)
        fixed.append(temp)
    return fixed

# Save regulons
def save_regulons(name):
    regulons_path = './regulons_' + name + '.csv'
    regs = pd.read_csv(regulons_path, header=2)
    regs = regs.rename(columns={'Unnamed: 8': 'TargetGenes'})
    reg_path = './' + name + '/'

    regs.drop(0, axis=0, inplace=True)
    regs.set_index('TF', inplace=True)

    a = create_target_genes(regs.TargetGenes.values)
    regs.loc[:, 'TargetGenes'] = a

    unique_idx = regs.index.unique()
    new_regs = {}
    for u in regs.index.unique():
        try:
            new_regs[u] = np.unique(
                np.concatenate([t for t in regs.loc[u, :].TargetGenes]))
        except:
            new_regs[u] = np.unique([t for t in regs.loc[u, :].TargetGenes])

    if not os.path.exists(reg_path):
        os.makedirs(reg_path)
    else:
        print('ERROR: Regulons already exist, ABORT!')
        return

    for u in regs.index.unique():
        with open(reg_path + u, 'w+') as f:
            for item in new_regs[u]:
                f.write("%s\n" % item)
