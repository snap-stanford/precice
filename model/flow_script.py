import pandas as pd
import numpy as np
import networkx as nx
import statsmodels.stats.api as sms
from scipy.stats import ncx2
from os.path import isfile
import os
from tqdm import tqdm
from flow import get_graph, get_expression_data,\
            add_weight, I_TF, get_TFs, solve,\
            solve_parallel, get_expression_lambda
import pdb

class flow_model():
    def __init__(self, args):
        self.args = args
        self.TFs = get_TFs(args.species)
        if args.expt_name == None:
            self.modelname = ('no_parallel_p_'+str(args.p_thresh) + 
                             '_max_step_'+str(args.max_step) +
                             '_hops_'+str(args.hops)+
                             '_DE_'+args.DE_data.split('/')[-1].split('.csv')[0])
            if args.regulon_name:
                self.modelname = self.modelname + '_regulon_' + args.regulon_name
            if args.adjacency:
                self.modelname = self.modelname + '_adjacency_' + args.adjacency.split('/')[-1].split('.csv')[0]
            if args.binary:
                self.modelname = 'binary_' + self.modelname
            if args.add_zero_TFs:
                self.modelname = 'inclZeroTFs_' + self.modelname
            if args.add_zero_genes:
                self.modelname = 'inclZeroGenes_' + self.modelname
            if args.pos_pert:
                self.modelname = 'PosPertOnly_' + self.modelname
            if args.pos_edges:
                self.modelname = 'PosEdges_' + self.modelname
            if args.top is not None:
                self.modelname = 'Top' +str(args.top) +'_' + self.modelname
            if len(args.remove_TFs) != 0:
                self.modelname = 'remove_TFs_' + '_'.join(self.args.remove_TFs)  +'_' + self.modelname

        else:
            self.modelname = args.expt_name
        print(self.modelname)

        # Set up graph structure
        if args.regulon_name is not None:
           G_df = get_graph(name = args.regulon_name,
                           TF_only=False, top=args.top)
           self.G = nx.from_pandas_edgelist(G_df, source=0,
                            target=1, create_using=nx.DiGraph())
        
        if args.adjacency is not None:
            self.read_weights = pd.read_csv(args.adjacency, index_col=0)
            self.G = nx.from_pandas_edgelist(self.read_weights, source='TF',
                            target='target', edge_attr = 'importance',
                            create_using=nx.DiGraph())
            G_df = self.read_weights

        # Remove TFs that are to be omitted
        self.TFs = list(set(self.TFs).difference(set(self.args.remove_TFs)))

        # Get differential expression data
        self.expression = get_expression_data(filename=args.DE_data,
                            graph_df=G_df, p_thresh=args.p_thresh,
                            zero_TFs=args.add_zero_TFs, TFs=self.TFs)

        try:
            self.expression = self.expression.loc[:,
                              ['gene_name', 'logFC']]  # limma
        except:
            self.expression = self.expression.loc[:,
                              ['gene_name', 'avg_log2FC']]  # Seurat

        # Get adjacency matrix
        self.adj_mat = self.create_adj_mat()

        self.B = self.expression.iloc[:, 1].values
        A = self.adj_mat.T
        if args.binary and args.pos_edges:
                A = np.array(A != 0).astype('float')

        # Set the diagonal elements to zero everywhere except the TFs
        if args.self_edges:
            print('Error: Not implemented self edges')
            return
        else:
            np.fill_diagonal(A, 0)
            each_hop = A.copy()
            last_hop = A.copy()
            for k in range(args.hops-1):
                last_hop = last_hop @ each_hop
                if args.binary:
                    A += last_hop/(k+2)
                else:
                    A += last_hop
        self.A = I_TF(A, self.expression, 1, self.TFs)

    def create_adj_mat(self):
        
        if self.args.regulon_name is not None:
           # For legacy code
           # Create a df version of the graph for merging
           G_df = pd.DataFrame(self.G.edges(), columns=['TF', 'target'])

           # Merge it with the weights DF
           self.read_weights = self.read_weights.merge(G_df, on=['TF', 'target'])
        
           for w in self.read_weights.iterrows():
              add_weight(self.G, w[1]['TF'], w[1]['target'], w[1]['importance'])

        if args.add_zero_genes:
            # Create a list of genes with zero differential expression
            TFs = self.TFs
            zero_DE_genes = [l for l in self.G.nodes()
                        if l not in self.expression['gene_name'].values]
            zero_DE_genes = [l for l in zero_DE_genes if l not in TFs]

            # Make the correction in the differential expression DF
            zero_DE_genes = pd.DataFrame(zero_DE_genes, columns=['gene_name'])
            colname = [c for c in self.expression.columns if c != 'gene_name'][0]
            zero_DE_genes[colname] = 0
            self.expression = self.expression.append(zero_DE_genes,
             ignore_index=True)

        # Get an adjacency matrix based on the gene ordering from the DE list
        return nx.linalg.graphmatrix.adjacency_matrix(
            self.G, nodelist=self.expression.gene_name.values).todense()

    def solve(self):

        # Create a list of lambdas
        lambdas = []

        if self.args.max_exp != False:
            for i in range(self.args.min_exp, self.args.max_exp):
                try:
                    lambdas += list(range(10 ** i, 10 ** (i + 1), 10 ** i))
                except:
                    lambdas = np.concatenate([lambdas,
                            np.arange(10 ** i, 10 ** (i + 1), 10 ** i)])
        elif self.args.max_step != 0:
            lambdas = np.arange(args.step_size, self.args.max_step,
                                args.step_size)
        else:
            print("Invalid input parameters to solver\n")
            return

        # Solve the system of equations
        if self.args.parallel:
            expression_lambda = solve_parallel(self.A, self.B, self.expression,
                                        lambdas, args.num_workers)
        else:
            expression_lambda = solve(self.A, self.B, self.expression, lambdas,
                                      positive=args.pos_pert)

        # Subset to TF's and save raw result
        TFs = self.TFs
        expression_lambda = expression_lambda[
            expression_lambda['gene_name'].isin(TFs)].reset_index(drop=True)
        expression_lambda.to_csv(os.path.join(self.args.out_dir, 
                          'perturbation_matrix_' + self.modelname + '.csv'))

    def k_active_sol(self, k, filtered=False, return_all_errors=False):
        # Given a perturbation matrix and a threshold,
        # find all solutions that result in k perturbations
        path = os.path.join(self.args.out_dir, 
                          'perturbation_matrix_' + self.modelname + '.csv')

        # Solve the system of linear equations if needed
        if isfile(path):
            pass
        else:
            self.solve()

        if k==0:
            return 'Source', self.get_error(np.zeros(len(self.B))), [' ']

        expression_lambda = get_expression_lambda(path)
        cols = [c for c in expression_lambda.columns
                                    if c not in ['gene_name', 'logFC']]
        TF_gene_names = expression_lambda.loc[:,'gene_name']
        expression_lambda = expression_lambda.loc[:, cols]
        total_genes_perturbed = (expression_lambda.abs() > 0).sum().values
        idx = np.where(total_genes_perturbed == k)[0]

        # Subsetting to lambda values that satisfy k condition
        if len(idx) == 0:
            return None, None, None
        k_perturb_df = expression_lambda.iloc[:, idx]

        k_perturb_dfs = []
        cols = [c for c in k_perturb_df.columns]
        k_perturb_df = k_perturb_df.assign(gene_name=TF_gene_names)

        for col in cols:
            k_perturb_dfs.append(k_perturb_df[k_perturb_df.loc[:, col].abs()
                            > 0].loc[:, ['gene_name', col]])
        if filtered:
            return self.filter_k_active_sol(k_perturb_dfs, return_all_errors)
        else:
            return k_perturb_dfs

    def get_error_interval(self, perturb_vec, sd=1):

        # Sampling method
        cov = np.diag((perturb_vec != 0).astype(int)) * (sd**2)
        samples = np.random.multivariate_normal(perturb_vec, cov, 1000)
        errors = [self.get_error(x) for x in samples]
        errors_iqr = [np.percentile(errors,25), np.percentile(errors,75)]
        return errors

        # TODO Closed form
        #chi2 = ncx2.fit(errors)
        #errors_ci = ncx2.interval(0.90, chi2[0], chi2[1], chi2[2], chi2[3])

        #mu = np.dot(self.A, perturb_vec)
        #sigma = np.dot(np.dot(self.A, cov), self.A.T)
        # The distribution over the perturbation outcomes can be obtained
        # This has to then be squared and summed
        # Estimate this as a non-central chi-square distribution


    def get_error(self, perturb_vec):
        error = np.sqrt(np.sum(np.square(self.B -
                                         np.matmul(self.A, perturb_vec))))
        return error


    def filter_k_active_sol(self, perturb_dfs, return_all_errors=False):
        # Measure the complete network perturbation produced
        # as well as its deviation from target
        # Return the perturbation set that causes minimum deviation

        perturb_errors = []
        opt_error = np.inf
        for perturb_df in perturb_dfs:
            perturb_vec = flow_model.create_perturbation_vec(
                self.expression, perturb_df).loc[:,'perturbation'].values
            error = self.get_error(perturb_vec)

            if error<opt_error:
                perturb_vec_save = perturb_vec
                opt_error = error
            # |y - (A*x)|
            perturb_errors.append(error)

        opt = np.argmin(perturb_errors)
        error_interval = self.get_error_interval(perturb_vec_save)
        if (return_all_errors):
            return perturb_dfs[opt], perturb_errors, error_interval

        else:
            return perturb_dfs[opt], perturb_errors[opt], error_interval

    def get_k_perturbation_range(self, max_range, save=False,
                                 return_all_errors=False):
        k_actives = {}
        print('Finding best solutions for each value of lambda:')
        for k in tqdm(range(max_range+1), desc='Processing'):
            sol, sol_error, error_interval = self.k_active_sol(k,
                            filtered=True, return_all_errors=return_all_errors)
            if isinstance(sol, str):
              if sol == 'Source':
                k_actives[k] = {'Genes': 'No perturbation',
                                'Perturbation': 0,
                                'Error':sol_error,
                                'Error reduction %': 0,
                                'Error_intervals':[]}
                source_error = float(sol_error)
            elif sol is not None:
                if len(sol) != 0:
                     k_actives[k] = {'Genes': sol.iloc[:, 0].values,
                                    'Perturbation': sol.iloc[:, 1].values,
                                    'Error':round(sol_error,3),
                                    'Error reduction %': 
                                        round((sol_error-source_error)/source_error,3),
                                    'Error_intervals':error_interval}
        result = pd.DataFrame(k_actives).T
        if save:
            result.to_pickle(os.path.join(self.args.out_dir, 
                          'active_' + self.modelname + '.pkl'))
        return result

    @staticmethod
    def create_perturbation_vec(express_df, perturb_df):
        output = express_df.loc[:, ['gene_name']]
        output = output.merge(perturb_df, how='outer', on='gene_name')
        output = output.fillna(0)
        cols = output.columns
        output = output.rename(columns={cols[0]: 'gene_name',
                                        cols[1]: 'perturbation'})
        return output


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Define solver arguments',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Parameters
    parser.add_argument('--expt_name', type=str,
                        default=None, help='(default: %(default)s)')
    parser.add_argument('--species', type=str, nargs='?',
                        default='human')

    # Lambda setting
    parser.add_argument('--max_step', type=int, default=51,
                        help='(default: %(default)s)')
    parser.add_argument('--step_size', type=float, default=1e-4,
                        help='(default: %(default)s)')
    parser.add_argument('--max_exp', type=int, default=False,
                        help='(default: %(default)s)')
    parser.add_argument('--min_exp', type=int, default=False,
                        help='(default: %(default)s)')
    parser.add_argument('--self_edges', type=bool, nargs ='?',
                        const=True, default=False)
    parser.add_argument('--hops', type=int, default=3,
                        help='(default: %(default)s)')

    # Solver parameters
    parser.add_argument('--parallel', type=bool, nargs ='?',
                        const=True, default=False)
    parser.add_argument('--num_workers', type=int,
                        default=80)
    parser.add_argument('--add_zero_genes', type=bool, nargs ='?',
                        const=True, default=False)
    parser.add_argument('--add_zero_TFs', type=bool, nargs ='?',
                        const=True, default=False)
    parser.add_argument('--p_thresh', type=float, default=1e-20,
                        help='(default: %(default)s)')
    parser.add_argument('--pos_pert', type=bool, nargs='?',
                        const=True, default=False, help='consider only positive'
                                                        ' perturbations')
    parser.add_argument('--binary', type=bool, nargs='?',
                        const=True, default=False)
    parser.add_argument('--pos_edges', type=bool, nargs='?',
                        const=True, default=False)
    parser.add_argument('--top', type=int, default=None,
                        help='Top percentile of edges to keep')
    parser.add_argument('--remove_TFs', nargs='+', default=[],
                        help='List of TFs to remove from consideration')

    # Input files and paths
    parser.add_argument('--adjacency', type=str,
                default=None)
    parser.add_argument('--regulon_name', type=str,
                        default=None)
    parser.add_argument('--DE_data', type=str,
                      default='DE_Friedman_1.csv',
                      help='(default: %(default)s)')
    parser.add_argument('--out_dir', type=str,
                        default='./')
    args = parser.parse_args()

    # Run flow script
    f = flow_model(args)

    # Identify k-active perturbations
    x = f.get_k_perturbation_range(10, save=True, return_all_errors=False)

    print('Done')
