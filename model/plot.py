import pandas as pd
import scanpy as sc
import numpy as np

import pandas as pd
import numpy as np
import seaborn as sns

import networkx as nx
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import matplotlib
matplotlib.rcParams['font.family'] = ['monospace']
matplotlib.rcParams['font.size'] = 18

from collections import defaultdict

from pathlib import Path

#Correcting plotting error where some rows do more than simply add one gene
def build_perturbation_labels(read_file):
    """
    For each row (after 'No perturbation' has been dropped),
    determine which gene(s) were added or dropped compared to the previous row.
    Returns a list of string labels like 'SNAI1 AND SOX2 AND NOT MESP2'.
    """
    labels = []
    prev_genes = set()
    import numpy as np

    for idx, row in read_file.iterrows():
        # Accept list, tuple, string, int...
        if isinstance(row.Genes, str):
            current_genes = set([row.Genes])
        elif isinstance(row.Genes, (list, tuple, np.ndarray)):
            # **Always cast every element to str**
            current_genes = set(str(g) for g in row.Genes)
        elif isinstance(row.Genes, (int, float)):
            # If just a single int/float (e.g., 0), treat as empty set
            current_genes = set()
        else:
            print("[DEBUG] Unknown Genes type in row", idx, type(row.Genes), row.Genes)
            current_genes = set()

        added = list(current_genes - prev_genes)
        dropped = list(prev_genes - current_genes)

        ### INFO BLOCK: print raw types for added/dropped ###
        print(f"[INFO] idx={idx}, row.Genes={row.Genes}, current_genes={current_genes}")
        print(f"[INFO]    added: {added} ({[type(a) for a in added]})")
        print(f"[INFO]    dropped: {dropped} ({[type(d) for d in dropped]})")

        # THE CRITICAL BIT: force cast all to str
        label_parts = []
        if added:
            label_parts.append(" AND ".join([str(x) for x in added]))
        if dropped:
            label_parts.append("NOT " + " AND NOT ".join([str(x) for x in dropped]))
        labels.append(" ".join(label_parts) if label_parts else "")
        prev_genes = current_genes
    return labels

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
    from pathlib import Path
    import os

    order = None
    if names is None:
        names = np.arange(len(active_files))

    active = {}
    pert_dist_dicts = {}
    for f, n in zip(active_files, names):
        p = Path(f)
        ext = p.suffix.lower()

        if ext == '.pkl':
            read_file = pd.read_pickle(p)
            if read_file['Genes'].values[0] == 'No perturbation':
                read_file = read_file[:max_rows]
                read_file = read_file.drop(0)
            active[n] = read_file

        elif ext == '.csv':
            # Use max_size if it's defined globally; otherwise call without it
            _max_size = globals().get('max_size', None)
            if _max_size is not None:
                active[n] = parse_csv(os.fspath(p), max_size=_max_size)
            else:
                active[n] = parse_csv(os.fspath(p))

        else:
            raise ValueError(f"Unsupported file extension: {ext}")

        if sort == False:
            order = get_rank_order(active[n])

        pert_dist_dicts[n] = get_pert_dist(active[n])
        active[f"{n}_labels"] = build_perturbation_labels(active[n])

    dicts = [pert_dist_dicts[n] for n in names]
    all_ticks = get_ticks(dicts)
    _ = [update_ticks(d, all_ticks) for d in dicts]

    print("[DEBUG]: dicts contents:")
    for i, d in enumerate(dicts):
        print(f"  dict {i}: {len(d)} genes")
        for k, v in d.items():
            print(f"    {k}: {v[:3]}{'...' if len(v) > 3 else ''}")

    labels_list = [active[f"{n}_labels"] for n in names]
    return dicts, labels_list, order

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

### INFO BLOCK to check if the input .pkl file has a "No perturbation" row, which it needs for downstream processing. This function called in make_precision_plot() (2025-06-03) ###
def has_no_perturbation_row(X):
    print("[INFO] has_no_perturbation_row called")  # <-- Indicates that the function is being called.
    print("[INFO] Columns in DataFrame:", X.columns)
    print("[INFO] First few rows:\n", X.head())

    likely_terms = ["none", "no perturbation", "no_perturbation", "null", "noperturbation"]
    found = False
    # Check all string/object columns for "no perturbation" (case-insensitive)
    for col in X.columns:
        if X[col].dtype == object or pd.api.types.is_string_dtype(X[col]):
            col_lower = X[col].astype(str).str.lower()
            for term in likely_terms:
                matches = col_lower.str.contains(term)
                if matches.any():
                    print(f'[INFO] Found a "no perturbation" row in column "{col}" (term: "{term}")')
                    print('[INFO] Rows:', X.loc[matches, :])
                    found = True
    if not found:
        print('[WARNING] No "no perturbation" row found (case-insensitive) in any string column!')
        # Optionally print unique values for inspection
        print('[INFO] First 10 unique values in first string column:', X.select_dtypes(include='object').iloc[:,0].unique()[:10])
    return found
### END INFO BLOCK ###  


#New version of make_precision_plot that solves erroneous arrow inheritance across iterations and erroenously showing only one arrow per composite label
def make_precision_plot(fname, max_rows: int = 12,
                        source: str = 'source', target: str = 'target',
                        out_dir: str | Path = "../data/plots",     # ← may be Path or str
                        out_name: str | None = None,               # ← file name only
                        ):
    """
    Plot cumulative precision and indicate, for each optimisation step,
    whether the newly added gene(s) should be up- (↑) or down-regulated (↓).

    Two fixes vs the original code:
      1. Arrow is based ONLY on the gene(s) introduced at this step,
         not on the average of the whole running list.
      2. If several genes are introduced in the same step, each gene in
         the label gets its own arrow (e.g.  "SOX2 ↑ AND MYC ↓").
    """
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import re

    # ----------------------------------------------------------
    #  Load file + helper-generated labels
    # ----------------------------------------------------------
    X = pd.read_pickle(fname)
    dicts, labels_list, order = create_pert_boxplot([fname],
                                                    sort=False,
                                                    max_rows=max_rows)
    labels = labels_list[0]          # only one active-file input
    precision_scores = (X['Error reduction %'][:max_rows] * -1).tolist()

    # ----------------------------------------------------------
    #  Decide arrows step-by-step
    # ----------------------------------------------------------
    perturbation_signs = [('None', '')]       # first row ⇒ “No perturb.”
    prev_genes: set[str] = set()

    for step_idx, label in enumerate(labels):
        row = X.iloc[step_idx + 1]            # skip “No perturbation”
        # normalise to list[str] / list[float]
        genes = ([row.Genes] if isinstance(row.Genes, str)
                 else list(row.Genes))
        pert = ([row.Perturbation] * len(genes)
                 if isinstance(row.Perturbation, (int, float))
                 else list(row.Perturbation))

        # which genes are NEW in this step?
        added = [g for g in genes if g not in prev_genes]
        dropped = [g for g in prev_genes if g not in genes]

        # ---------- build label ------------------------------------------------
        parts_added = [
            f"{g}{'↑' if pert[genes.index(g)] > 0 else '↓'}"
            for g in added
        ]

        label_parts = parts_added.copy()               # start with added genes
        if dropped:                                    # tack on any removals
            label_parts.append(
                "NOT " + " AND NOT ".join(map(str, dropped))
            )

        label_with_arrows = " AND ".join(label_parts) if label_parts else label

        perturbation_signs.append((label_with_arrows, ""))
        prev_genes = set(genes)        # refresh current-gene set

    # keep both lists the same length
    precision_scores = precision_scores[:len(perturbation_signs)]

    # ----------------------------------------------------------
    #  Plotting
    # ----------------------------------------------------------
    import pandas as pd
    df = pd.DataFrame({'Precision Score': precision_scores,
                       'Perturbation': perturbation_signs})

    plt.figure(figsize=(5, 8))
    pos = range(len(df))
    plt.plot(df['Precision Score'], pos,
             color='purple', linewidth=4, zorder=1)
    plt.scatter(df['Precision Score'], pos,
                color='white', edgecolor='black',
                s=400, zorder=2)

    ax = plt.gca()
    ax.set_xlabel('Precision Score (%)', fontweight='bold')
    ax.set_title(f'{source} → {target}')

    y_labels = []
    for i, (lbl, _) in enumerate(df['Perturbation']):
        if i >= 2:
            y_labels.append("AND " + lbl)
        else:
            y_labels.append(lbl)
    ax.set_yticks(pos)
    ax.set_yticklabels(y_labels, ha='right')

    ax.invert_yaxis()
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xticklabels(['0', '10', '20', '30', '40', '50'])

    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(2)
        ax.spines[spine].set_color('black')

    plt.grid(True, axis='y', linestyle='--', color='grey')

    ### New: allow specifying file directory and -name (2025-07-09):
    from pathlib import Path
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # if caller gives a name, use it; otherwise fall back to the old pattern
    if out_name is None:
        out_name = f"{source}_to_{target}.jpeg"

    out_path = out_dir / out_name
    ###
    #Old lines that don't allow specifying file name or path (commented out 2025-07-09)
    #out_dir = "../data/plots"
    #os.makedirs(out_dir, exist_ok=True)
    #out_path = f"{out_dir}/{source}_to_{target}.jpeg"
    
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.show()

    return out_path