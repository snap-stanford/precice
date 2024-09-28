import sys
import pandas as pd
import numpy as np
import networkx as nx

sys.path.append('../model/')
from flow import get_graph
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.neural_network import MLPRegressor

no_model_count = 0

def nonzero_idx(mat):
    mat=pd.DataFrame(mat)
    return mat[(mat > 0).sum(1) > 0].index.values

def data_split(X, y):
    nnz = list(set(nonzero_idx(X)).intersection(set(nonzero_idx(y))))

    if len(nnz) <= 1:
        global no_model_count
        no_model_count += 1

        return -1,-1

    train_split, val_split = train_test_split(nnz, test_size=0.15)
    return train_split, val_split

def train_regressor(X, y, kind='linear'):

    if kind == 'linear':
        model = linear_model.LinearRegression()
    elif kind == 'lasso':
        model = linear_model.Lasso(alpha=10e4)
    elif kind == 'elasticnet':
        model = linear_model.ElasticNet(alpha=10e4, l1_ratio=0.5,
                                        max_iter=1000)
    elif kind == 'ridge':
        model = linear_model.Ridge(alpha=10e4, max_iter=1000)
    elif kind == 'MLP':
        model = MLPRegressor(hidden_layer_sizes=(20,10,5), max_iter=1000)

    reg = model.fit(X, y)
    loss = np.mean(y - model.predict(X.values))
    return reg, loss, reg.score(X, y)


def evaluate_regressor(model, X, y):
    y_cap = model.predict(X.values)
    loss = np.mean(y - y_cap)

    return loss

def init_dict():
    d = {}
    d['linear'] = []
    d['ones'] = []
    d['lasso'] = []
    d['elasticnet'] = []
    d['ridge'] = []
    d['MLP'] = []
    return d

def get_weights(adj_mat, X, lim=20000):
    models = init_dict()
    val_loss = init_dict()
    train_loss = init_dict()
    train_score = init_dict()

    adj_mat_idx = np.arange(len(adj_mat))
    new_adj_mat = np.zeros(adj_mat.shape)
    np.random.shuffle(adj_mat_idx)
    count = 0

    def trainer(kind, feats, y, train_split, val_split):
        model, train_loss_, train_score_ = train_regressor(
                                        feats.loc[train_split,:],
                                        y.loc[train_split], kind=kind)
        val_loss_ = evaluate_regressor(model,
                                       feats.loc[val_split, :],
                                       y.loc[val_split])

        # Store results
        val_loss[kind].append(val_loss_)
        train_loss[kind].append(train_loss_)
        train_score[kind].append(train_score_)
        try: models[kind].append(model.coef_);
        except: pass;

    for itr in adj_mat_idx:
        i = adj_mat[itr]
        if i.sum() > 0:
            idx = np.where(i > 0)[1]

            feats = X.iloc[:, idx]
            y = X.iloc[:, itr]
            train_split, val_split = data_split(feats, y)

            if train_split==-1:
                continue

            print(count)

            # Linear Regression models
            trainer('linear', feats, y, train_split, val_split)
            #trainer('ridge', feats, y, train_split, val_split)
            #trainer('MLP', feats, y, train_split, val_split)
            #trainer('elasticnet', feats, y, train_split, val_split)

            # All edges are 1
            model = linear_model.LinearRegression()
            model.coef_ = np.ones(len(idx))
            model.intercept_ = 0
            val_loss_ = evaluate_regressor(model, feats.loc[val_split,:],
                                                  y.loc[val_split])
            val_loss['ones'].append(val_loss_)

            # Add row to new weight matrix
            for j,k in enumerate(idx):
                new_adj_mat[itr][k] = models['linear'][-1][j]

            count += 1

        if count >= lim:
            break

    return train_loss, train_score, val_loss, models, new_adj_mat

regulons_path = '../Data/regulons_MEF_epicardial.csv'
MEF_epicardial = pd.read_csv('../Data/MEF_epicardial.txt', delimiter='\t')
regs = pd.read_csv(regulons_path, header=2)

G = get_graph(name = 'MEF_epicardial',
                   TF_only=False)
G = nx.from_pandas_edgelist(G, source=0,
                    target=1, create_using=nx.DiGraph())
adj_mat = nx.linalg.graphmatrix.adjacency_matrix(G).todense().T

# Remove self-edges
np.fill_diagonal(adj_mat, 0)

X = MEF_epicardial.iloc[:,1:].T

# Train and test index split
train_idx, test_idx = train_test_split(X.index.values, test_size=0.10)
train_data = X.loc[train_idx, :]

train_loss, train_score, val_loss, models, new_adj_mat = get_weights(adj_mat,
                                            train_data, lim=20000)

# Save final results
np.save('train_loss', train_loss)
np.save('val_loss', train_loss)
np.save('linear_adj_mat', new_adj_mat)

# Convert coefficients into new weight matrix
print('Done')
