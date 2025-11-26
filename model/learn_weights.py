import sys
import pandas as pd
import numpy as np
import scanpy as sc
import networkx as nx
sys.path.append('../model/')
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.neural_network import MLPRegressor
import argparse

no_model_count = 0

def nonzero_idx(mat):
    mat=pd.DataFrame(mat)
    return mat[(mat > 0).sum(1) > 0].index.values

def data_split(X, y, size=0.1):
    nnz = list(set(nonzero_idx(X)).intersection(set(nonzero_idx(y))))

    if len(nnz) <= 2:
        global no_model_count
        no_model_count += 1

        return -1,-1

    train_split, val_split = train_test_split(nnz, test_size=size,
                                              random_state=42)
    return train_split, val_split

def train_regressor(X, y, kind, alpha=10): #Default alpha=0.

    if kind == 'linear':
        model = linear_model.LinearRegression()
    elif kind == 'lasso':
        model = linear_model.Lasso(alpha=alpha)
    elif kind == 'elasticnet':
        model = linear_model.ElasticNet(alpha=alpha, l1_ratio=0.5,
                                        max_iter=1000) #Default max_iter=1000. 
    elif kind == 'ridge':
        model = linear_model.Ridge(alpha=alpha, max_iter=1000) #Default max_iter=1000 
    elif kind == 'MLP':
        model = MLPRegressor(hidden_layer_sizes=(20,10), max_iter=1000) #Default max_iter=1000 

    reg = model.fit(X, y)
    loss = np.sqrt(np.mean((y - model.predict(X))**2))
    return reg, loss, reg.score(X, y)


def evaluate_regressor(model, X, y):
    y_cap = model.predict(X)
    loss = np.sqrt(np.mean((y - y_cap)**2))

    return loss, y, y_cap

def init_dict():
    d = {}
    d['linear'] = []
    d['lasso'] = []
    d['ridge'] = []
    d['MLP'] = []
    return d

# Looks at the median of max expression across cells/not genes
def max_median_norm(df):
    return df/df.max().median()

def get_weights(adj_mat, exp_adata, nodelist, method='linear', lim=50000):
    models = init_dict()
    adj_list = {}

    adj_list['TF'] = []; adj_list['target'] = []; adj_list['importance'] = [];

    adj_mat_idx = np.arange(len(adj_mat))
    np.random.shuffle(adj_mat_idx)
    count = 0


    def trainer(kind, feats, y):
        model, _, _ = train_regressor(
                                        feats, y, kind=kind)

        # Store results
        try: models[kind].append(model.coef_);
        except: pass;


    def trainer_split(kind, feats, y, train_split, val_split):
        models_ = []
        val_losses_ = []
        
        for alpha in [1e-6, 1e-4, 1e-2, 1e-1]:
             model, _, _ = train_regressor(
                                        feats[train_split,:],
                                        y[train_split], kind=kind,
                                        alpha=alpha)
             val_loss, _, _ = evaluate_regressor(model,
                                       feats[val_split, :],
                                       y[val_split])

             models_.append(model)
             val_losses_.append(val_loss)
        
        best_model = models_[np.argmin(val_losses_)]

        # Store results
        try: models[kind].append(best_model.coef_);
        except: pass;


    print('T genes: ', str(len(adj_mat_idx)))
    for itr in adj_mat_idx:
        i = adj_mat[itr]
        #print("INFO: itr =", itr)
        #print("INFO: i shape:", i.shape)
        #print("INFO: i contents:", i)
        # If sum is zero, print and continue
        if i.sum() > 0:
            ### INFO BLOCK ###
            #print("INFO: Shape of i:", i.shape)
            #print("INFO: Result of np.where(i > 0):", np.where(i > 0))
            ### END INFO BLOCK ###
            idx = np.where(i > 0)[0] #Default: idx = np.where(i > 0)[1] but i is a 1D array (the input csv file is a 2D matrix but in it, each row (representing an edge) is a 1D array, so cannot be subsetted for [1]
            TFs = np.array(nodelist)[idx]
            target = np.array(nodelist)[itr]

            feats = exp_adata[:, TFs].X.toarray()
            y = exp_adata[:, target].X.toarray()
            
            if method=='linear': 
                trainer('linear', feats, y)
            else:
                train_split, val_split = data_split(feats, y, size=0.1)
                if train_split==-1: continue;
                trainer_split(method, feats, y, train_split, val_split)            

            # Add row to new weight matrix
            for j,k in enumerate(TFs):
                adj_list['TF'].append(k)
                adj_list['target'].append(target)
                try:
                    adj_list['importance'].append(models[method][-1][0][j])
                except:
                    adj_list['importance'].append(models[method][-1][j])

            print(count)
            count += 1

        if count >= lim:
            break
    return models, adj_list


### HELPER FUNCTION THAT ATTEMPTS TO HARMONIZE GENE NAMES BETWEEN THE [...]LINEARWEIGHTS.CSV FILE AND THE ANNDATA FILE (e.g., CXORF40A-type names, which are CXorf40A in the AnnData) 2025-06-25###
def harmonize_gene_names(df, adata, cols=['TF', 'target']):
    """
    Attempt to harmonize gene names in the network CSV with those in AnnData.
    For each column (TF, target), if a gene is not found in adata.var_names,
    look for a case-insensitive match (e.g., C5ORF24 → C5orf24).
    """
    import difflib

    varnames = np.array(adata.var_names)
    varnames_lower = {name.lower(): name for name in varnames}

    for col in cols:
        new_names = []
        for g in df[col]:
            # Try exact match
            if g in varnames:
                new_names.append(g)
            # Try case-insensitive match
            elif g.lower() in varnames_lower:
                new_names.append(varnames_lower[g.lower()])
            else:
                # Try fuzzy match (optional)
                matches = difflib.get_close_matches(g, varnames, n=1, cutoff=0.8)
                if matches:
                    print(f"[WARN] '{g}' replaced by closest match '{matches[0]}'")
                    new_names.append(matches[0])
                else:
                    print(f"[ERROR] Could not match gene '{g}' in AnnData. It will be kept as is (may cause error later).")
                    new_names.append(g)
        df[col] = new_names
    return df

def main(args):
    from pathlib import Path

    # Coerce inputs to pathlib.Path
    data_path  = Path(args.data_path)
    graph_path = Path(args.graph_name)
    out_dir    = Path(args.out_dir)

    # Read inputs
    adata = sc.read_h5ad(str(data_path))
    G = pd.read_csv(str(graph_path), header=0)  # keep header=0 per your note

    # Harmonize gene names
    G = harmonize_gene_names(G, adata, cols=['TF', 'target'])

    # Build directed graph and adjacency
    G_nx = nx.from_pandas_edgelist(G, source='TF', target='target', create_using=nx.DiGraph())
    adj_mat = nx.linalg.graphmatrix.adjacency_matrix(G_nx).todense().T
    nodelist = [n for n in G_nx.nodes()]

    # Remove self-edges
    np.fill_diagonal(adj_mat, 0)

    # Learn weights
    models, adj_list = get_weights(adj_mat, adata, nodelist, method=args.method, lim=20000)

    # Prepare output path and write
    out_dir.mkdir(parents=True, exist_ok=True)
    out_name = graph_path.stem
    out_file = out_dir / f"{out_name}_{args.method}_learntweights.csv"
    pd.DataFrame(adj_list).to_csv(out_file, index=False)

    print(f"Done. Saved learnt weights to: {out_file}")
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Set model hyperparametrs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #torch.cuda.set_device(4)

    parser.add_argument('--data_path', type=str,
                        help='Directory for AnnData')
    parser.add_argument('--graph_name', type=str,
                        help='Graph filename')
    parser.add_argument('--out_dir', type=str,
                        help='Output filename')
    parser.add_argument('--method', type=str,
                        help='Regression method')


    parser.set_defaults(
    data_path ='../Notebooks/Friedman_1.h5ad',
    graph_name='../Data/transcription_networks/G_all_edges_Friedman_1',
    method='linear',
    out_dir='./')

    args = parser.parse_args()
    main(args)
