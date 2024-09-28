import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import kstest, ks_2samp


# Given a TF, find all edges linked to its regulon
def find_regulon_edges(TF, mappings, verbose=True):
    try:
        regulon_members = pd.read_csv('../Data/regulons/' + TF + '.txt', header=None).iloc[:, 0].values
    except:
        if (verbose):
            print(TF + ':File does not exist')
        return
    TF_id = mappings[TF]
    regulons_id = []
    for r in regulon_members:
        try:
            regulons_id.append(mappings[r])
        except:
            pass
    return [(TF_id, r) for r in regulons_id]


# Assign edge weight to a single edge
def add_weight(G, u, v, weight):
    try:
        G.remove_edge(u, v)
    except:
        # If the edge doesn't exist don't add a weighted version
        # return
        pass
    G.add_edge(u, v, weight=weight)


def get_max_out_deg(g, weight=False):
    if (weight):
        return list(g.nodes)[np.argmax([g.out_degree(i, weight='weight') for i in g.nodes])]
    else:
        return list(g.nodes)[np.argmax([g.out_degree(i) for i in g.nodes])]


# Given edge list, add corresponding weights to edges
def add_weight_regulon(G, edges, weight):
    for u, v in edges:
        add_weight(G, u, v, weight)


# Using shannon diversity identify and rank nodes needed to break down network
# Given a graph, calculate its normalized Shannon diversity
def shannon(g):
    p_i = [len(x) / len(g) for x in nx.weakly_connected_components(g)]
    return -np.sum([(p * np.log(p)) for p in p_i]) / np.log(len(g))


# Identify the node the loss of which maximizes increase in Shannon diversity
def shannon_max_step(g, steps=None):
    ents = []

    if steps == None:
        steps = g.nodes()
    for n in steps:
        g_ = g.copy()
        g_.remove_node(n)
        ents.append([n, shannon(g_)])
    return np.sort(ents, 0)[-1][0]


# Return ranked nodes along with their induced entropies
# and number of remaining wccs
def get_ranked_list(g, method='shannon', max_itr=200):
    g_ = g.copy()
    ents = [shannon(g)]
    wccs = [nx.number_weakly_connected_components(g)]
    nodes = [None]

    for it in range(len(g_)):
        if method == 'shannon':
            drop = shannon_max_step(g_)
        elif method == 'max_out_deg':
            drop = get_max_out_deg(g_, weight=True)
        print(it)
        g_.remove_node(drop)

        nodes.append(drop)
        ents.append(shannon(g_))
        wccs.append(nx.number_weakly_connected_components(g_))

        if (it==max_itr):
            break

    return (nodes, ents, wccs)


# Given AUC vector, add corresponding weights to existing graph
def AUC_to_weighted_graph(G, vec, mappings):
    G1 = nx.Graph.copy(G)
    for TF, weight in vec.iteritems():
        TF = TF.split('(')[0]
        edges = find_regulon_edges(TF, mappings)
        if edges is not None:
            add_weight_regulon(G1, edges, weight)
    return G1


class AUCell:
    def __init__(self, cell_type, df, type_dict):
        self.df = df

        if type(cell_type) is not list:
            self.cell_type = [cell_type]
        else:
            self.cell_type = cell_type
        self.IDs = np.concatenate([type_dict[l] for l in self.cell_type])

    # Returns the p-value for the KS test between two distributions
    @staticmethod
    def d_AUCell(AUCell_1, AUCell_2, gene, df=None):
        if df is None:
            df = AUCell_1.df
        x = df.loc[AUCell_1.good_idx(AUCell_1.IDs), gene].values
        y = df.loc[AUCell_2.good_idx(AUCell_2.IDs), gene].values
        return ks_2samp(x, y)[1]

    # Run the KS test across all TF for given cell types
    def d_AUCell_allTF(self):
        self.data = {gene: self.d_AUCell(gene) for gene in self.df.columns}

    # Use only indices that are in the dataframe
    def good_idx(self, list_):
        return [s for s in list_ if s in self.df.index]
