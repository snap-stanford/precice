import pandas as pd
import scanpy as sc
import numpy as np

import pandas as pd
import numpy as np
import seaborn as sns

import networkx as nx

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import matplotlib
font = {'size'   : 18}
matplotlib.rc('font', **font)
from collections import defaultdict

def get_pert_sums(df):
    gene_pert = defaultdict(int)

    for k,v in df.iterrows():
        for it, g in enumerate(v['Genes']):
            gene_pert[g] += v['Perturbation'][it]
            
    return gene_pert

def get_pert_dist(df):
    gene_pert = {}
    
    for k,v in df.iterrows():
        for it, g in enumerate(v['Genes']):
            try:
                gene_pert[g].append(v['Perturbation'][it])
            except:
                if type(v['Perturbation']) == int:
                    continue
                else:
                    gene_pert[g] = [v['Perturbation'][it]]
            
    return gene_pert

def get_rank_order(in_df):
    current_set = set()
    order = []
    for x in in_df['Genes']:
        s = set(x)
        order.extend([item for item in set(x).difference(current_set)])
        current_set = current_set.union(set(x))
        
    return order[::-1]

def plot_overall_gene_sums(dict_):
    plt.barh(list(dict_.keys()), list(dict_.values()), alpha=1)
    plt.plot([0,0],[-1,len(dict_)+1], 'k')
    
def plot_overall_gene_box_plot(dict_, color='red', order=None, ax=None):
    plot_df = pd.DataFrame.from_dict(dict_, orient='index')
    if order is None:
        plot_df = plot_df.sort_index()
    else:
        plot_df = plot_df.loc[order,:]
    plot_df = plot_df.melt(ignore_index=False).dropna()
    plot_df = plot_df.reset_index()

    plot_vals = [plot_df[plot_df['index'] == p]['value']
                     for p in plot_df['index'].unique()]
    
    if ax == None:
        #ax = sns.boxplot(data = plot_df, x='value', y='index', color=color)
        bp = plt.boxplot(plot_vals, vert=False, showfliers=False, 
                         showmeans=False, patch_artist=True, widths=0.65, alpha=1)
        ax = plt.gca()
    else:
        #sns.boxplot(data = plot_df, x='value', y='index', color=color, ax=ax)
        bp = plt.boxplot(plot_vals, vert=False, showfliers=False, 
                         showmeans=False, patch_artist=True, widths=0.65)
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color='black', alpha=1)
    
    for patch in bp['boxes']:
        patch.set(facecolor=color) 
        
    plt.plot([0,0],[-1,len(dict_)+1], 'k')
    plt.gca().grid(axis='y', linestyle='--')
    plt.ylabel('Predicted TF to perturb')
    plt.xlabel('Magnitude of perturbation')
    
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .8))
    
    return ax
    
def get_ticks(dicts):
    all_keys = set()
    for d in dicts:
        all_keys = all_keys.union(set(d.keys()))
    return all_keys

def update_ticks(dict_, all_ticks):
    for k in all_ticks:
        if k not in dict_.keys():
            dict_[k] = [0]
            
def create_pert_boxplot(active_files, names=None, sort=True, max_rows=12):
    
    order = None
    if names is None:
        names = np.arange(len(active_files))
    
    active = {}
    pert_dist_dicts = {}
    for f,n in zip(active_files, names):
        if f[-3:] == 'pkl':
            read_file = pd.read_pickle(f)
            if read_file['Genes'].values[0] == 'No perturbation':
                read_file = read_file.drop(0)
                read_file = read_file[:max_rows]
            active[n] = read_file
        elif f[-3:] == 'csv':
            active[n] = parse_csv(f, max_size=max_size)
        if sort == False:
            order = get_rank_order(active[n])
        pert_dist_dicts[n] = get_pert_dist(active[n])
        
    dicts = [pert_dist_dicts[n] for n in names]
    all_ticks = get_ticks(dicts)
    _ = [update_ticks(d, all_ticks) for d in dicts]
    
    return dicts, order

def plot_pert_boxplot(active_files, names=None, dicts=None, title='None',
                     sort=True, color=None, ax=None):
    if ax == None:
        plt.figure(figsize = [10,10])

    if color is None:
        colors = ['red',
             'green',
             'blue',
             'yellow',
             'purple']
    else:
        colors = [color]
    axs = []
    
    if names is None:
        names = np.arange(len(active_files))
        order = None
        
    if dicts is None:
        dicts, order = create_pert_boxplot(active_files, names, sort=sort)
        
    for i,d in enumerate(dicts):
        if ax == None:
            axs.append(plot_overall_gene_box_plot(d, colors[i], 
                                                  order=order))
        else:
            plot_overall_gene_box_plot(d, colors[i], order=order, ax=ax)
            axs = [ax]

    plt.title(title)

    #Make custom legend
    legend_elements = [Patch(facecolor=c, 
                             edgecolor=c,
                             alpha = 1,
                             label=d,
                            ) for c,d in zip(colors[:len(names)], names)]

    # Create the figure
    axs[-1].legend(handles=legend_elements) 
    
    # Is this guaranteed to be in the right order
    return order

def make_bar_plot(fname, name):
    fnames = [fname]

    names = [name]
    #names = ['Endoderm'] 
             #'Mesoderm', 'Endoderm+Mesoderm']

    fig = plt.figure(figsize=[10,10])
    ax1 = plt.gca()

    gene_list = plot_pert_boxplot(fnames, names, sort=False, title='',
                      color='forestgreen', ax=ax1)
    ax1.set_ylim([0,len(gene_list)+1])
    ax1.set_xlim([-8,8])

    ax1.set_yticklabels(gene_list)

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.spines[axis].set_color("black")

    #sns.set_style("whitegrid", {'grid.linestyle': '--'})
    #plt.gca().grid(axis='x')

    plt.xlabel(r'Magnitude of perturbation, $\theta$')

    plt.ylabel(r'$\leftarrow$ Ranked TFs to perturb')
    
    return gene_list
    
    
def make_error_plot(fname, max_rows=12):
    
    X = pd.read_pickle(fname)
    plt.figure(figsize = [3,10])
    pos = [len(g) for g in X['Genes'][1:max_rows]]
    plt.plot(X['Error reduction %'][:max_rows], [0] + pos, 
             color='forestgreen', linewidth=2)
    plt.scatter(X['Error reduction %'][:max_rows], [0] + pos,  
             color='forestgreen')
    ax = plt.gca()
    plt.xlabel('% Error reduction \n over no perturb. \n (cumulative)')

    ax.set_xticks([0, -0.1,  -0.2, -0.3, -0.4, -0.5 ])
    ax.set_xticklabels(['0','-10','-20', '-30', '-40', '-50'])

    ax.set_xlim([-0.65,0.05])
    ax.set_ylim([0,max_rows])
    plt.yticks(range(max_rows))
    ax.set_yticklabels('')

    #sns.set_style("whitegrid", {'grid.linestyle': '--'})
    plt.gca().grid(axis='x')

    plt.gca().invert_yaxis()

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
        ax.spines[axis].set_color("black")
