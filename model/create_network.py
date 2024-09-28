import pandas as pd
import numpy as np
import glob as glob
import networkx as nx
import utils
from flow import get_TFs

def create_net(args):
    if args.name is None:
        args.name = args.regulons.split('/')[-1]

    # First convert edge lists into adjacency matrices
    edges = []
    for r in glob.glob(args.regulons + '/*'):
        read = pd.read_csv(r, header=0, names=[1])
        read[0] = r.split('/')[-1].split('.')[0]
        edges.append(read)

    # Put all the list together and
    # create a mapping of every unique gene
    edges = pd.concat(edges).reset_index(drop=True)
    all_genes = np.unique(edges.values.flatten())
    mappings = {all_genes[i]: i for i in range(len(all_genes))}

    # Re-read edges with idx values instead of gene names
    # This is repetitive but works faster than editing the df
    edges = []
    for r in glob.glob(args.regulons + '/*'):
        read = pd.read_csv(r, header=0, names=[1])
        read[1] = [mappings[k] for k in read[1].values]
        read[0] = mappings[r.split('/')[-1].split('.')[0]]
        edges.append(read)

    # Put all the list together
    edges = pd.concat(edges).reset_index(drop=True)
    edges = edges.loc[:,[0,1]]
    edges = edges.rename(columns={0:'source',1:'target'})

    # Save as graph
    G = nx.from_pandas_edgelist(edges,create_using=nx.DiGraph())
    G_relabelled = nx.relabel_nodes(G,utils.invert_dict(mappings))
    nx.write_edgelist(G_relabelled, 'G_all_edges_'+args.name,
                      delimiter=',', data=dict(G.out_degree(G.nodes())))

    # Remove all non-TF genes
    TF_nodes = get_TFs(args.species)
    mappings_T = utils.invert_dict(mappings)
    node_names = list(G.nodes())
    _ = [G.remove_node(i) for i in node_names if mappings_T[i] not in TF_nodes]

    # Save new graph
    G_relabelled = nx.relabel_nodes(G,utils.invert_dict(mappings))
    nx.write_edgelist(G_relabelled, 'TF_only_'+args.name,
                      delimiter=',', data=dict(G.out_degree(G.nodes())))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Define solver arguments',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Input files
    parser.add_argument('--regulons', type=str, default='../Data/regulons_MEF')
    parser.add_argument('--species', type=str, default='human')
    parser.add_argument('--name', type=str, default=None)
    args = parser.parse_args()

    create_net(args)
