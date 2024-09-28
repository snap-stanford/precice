## Script for visualizing networks

import graphviz
import pygraphviz
from networkx.drawing.nx_agraph import to_agraph

import colorsys
import matplotlib as mpl


def float_to_hsv(val, norm, cmap):
    m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mm = m.to_rgba(val)
    return colorsys.rgb_to_hsv(mm[0], mm[1], mm[2])


def get_TF_expression_data():
    # Read in normalized expression data
    expression_data = pd.read_csv('../Data/normed_data_all_ss2_hidden_names.txt', delimiter='\t')

    # Read in TF network
    df = pd.read_csv('../Data/TF_only', header=None)

    expression_genes = np.unique(np.array(df.values).flatten())
    TF_expression = expression_data[expression_data['gene_name'].isin(expression_genes)]

    return TF_expression


## Create a cell type specific sub-df for plotting
def get_plot_df(in_df, cell_type, func='mean'):
    plot_df = in_df[in_df['variable'] == cell_type].loc[:, ['gene_name', 'value']]
    plot_df['value'] = plot_df.value.astype('float')
    if func == 'mean':
        plot_df = plot_df.groupby('gene_name').mean()
    elif func == 'median':
        plot_df = plot_df.groupby('gene_name').median()
    return plot_df


def get_graph():
    df = pd.read_csv('../Data/TF_only', header=None)
    G = nx.from_pandas_edgelist(df, source=0, target=1, create_using=nx.DiGraph())
    return G


def get_AGraph(G):
    A = to_agraph(G)

    # set some default node attributes
    A.node_attr['style'] = 'filled'
    A.node_attr['shape'] = 'circle'
    A.node_attr['fixedsize'] = 'false'
    A.node_attr['fontcolor'] = 'white'
    A.node_attr['color'] = 'black'
    A.graph_attr['overlap'] = 'false'

    return A


def draw_graphviz(G, cell_type, norm, df_):
    cmap = mpl.cm.Oranges
    for it, i in enumerate(G.nodes()):
        # express = int(np.max([1,plot_df.loc[i].value]))
        M = float_to_hsv(np.max([0, df_.loc[i].value]), norm, cmap)
        if df_.loc[i].value == 0:
            A.get_node(i).attr['fillcolor'] = "#ffffff"
            A.get_node(i).attr['fontcolor'] = 'black'
        else:
            A.get_node(i).attr['fillcolor'] = "%f, %f, %f" % (M[0], M[1], M[2])
        # A.get_node(i).attr['fillcolor']="#%2x0000"%(express*1)

    A.layout('dot')
    # A.graph_attr.update(fontsize=10)
    A.draw(cell_type)


TF_expression = get_TF_expression_data()
bars = TF_expression.melt(id_vars='gene_name')
bars = bars[bars['variable'] != 'gene_id']
bars['variable'] = [v.split('_')[0] for v in bars['variable']]

G = get_graph()
A = get_AGraph(G)

cell_types = ['inhib', 'excite', 'stem']

plot_dfs = []
for c in cell_types:
    df_ = get_plot_df(bars, c)
    plot_dfs.append(df_.apply(lambda x: np.log(x+1e-5)))

all_plot_dfs = pd.concat(plot_dfs)
norm = mpl.colors.Normalize(vmin = np.min(all_plot_dfs.values), vmax = np.max(all_plot_dfs.values))

for it,c in enumerate(cell_types):
    draw_graphviz(G, c, norm, plot_dfs[it])

## Box plot
order = stem_df.groupby('gene_name').mean().sort_values('value', ascending=False).index.values

fig, ax = plt.subplots(figsize=[6,50])
plt.xscale('log')
sns.boxplot(y='gene_name', x='value', hue='variable', data=bars, order=order)