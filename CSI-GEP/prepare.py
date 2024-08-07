import pandas as pd
from os.path import join as pjoin
import os, errno, gc
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'
import argparse
import numpy as np
import scanpy as sc
from scipy import io, sparse
from random import choices

def counts_split(raw_counts_path, output_path, numgenes=2000):
    """
    Split input raw count matrix into 2 halves and write to h5ad files in the current directory
    Compute tpm and save to output folder
    Find HVGs and write to output folder
    Normalize count data and save to output folder
    """
    check_dir_exists(output_path)
    check_dir_exists(pjoin(output_path, "cNMF_split0"))
    check_dir_exists(pjoin(output_path, "cNMF_split1"))
    check_dir_exists(pjoin(output_path, 'cNMF_split0/cnmf_tmp'))
    check_dir_exists(pjoin(output_path, 'cNMF_split1/cnmf_tmp'))

    data = sc.read_h5ad(raw_counts_path)
    splits = np.array(choices(["split0", "split1"], k=data.n_obs))
    adata_split0 = data[splits=="split0"]
    adata_split1 = data[splits=="split1"]
    del adata_split0.raw
    del adata_split1.raw

    adata_split0.write_h5ad(pjoin(output_path, "split0.edited.h5ad"),compression='gzip')
    adata_split1.write_h5ad(pjoin(output_path, "split1.edited.h5ad"),compression='gzip')
    
    # split0
    print("Computing TPM values on split0...")
    tpm = compute_tpm(adata_split0)
    sc.write(pjoin(output_path, 'cNMF_split0/cnmf_tmp/cNMF_split0.tpm.h5ad'), tpm)
    gene_tpm_mean = np.array(tpm.X.mean(axis=0)).reshape(-1)
    gene_tpm_stddev = var_sparse_matrix(tpm.X)**.5
    input_tpm_stats = pd.DataFrame([gene_tpm_mean, gene_tpm_stddev], index = ['_mean', '_std']).T
    save_df_to_npz(input_tpm_stats, pjoin(output_path, 'cNMF_split0/cnmf_tmp/cNMF_split0.tpm_stats.df.npz'))

    print("Computing highly variable genes and normalized counts on split0...")
    highvargenes = None
    genels_directory = pjoin(output_path, 'cNMF_split0/cNMF_split0.overdispersed_genes.txt')

    norm_counts = get_norm_counts(adata_split0, tpm, genels_directory, 
     high_variance_genes_filter=highvargenes, num_highvar_genes=numgenes)

    sc.write(pjoin(output_path, 'cNMF_split0/cnmf_tmp/cNMF_split0.norm_counts.h5ad'),norm_counts)
    
    # split1
    print("Computing TPM values on split1...")
    tpm = compute_tpm(adata_split1)
    sc.write(pjoin(output_path, 'cNMF_split1/cnmf_tmp/cNMF_split1.tpm.h5ad'), tpm)
    gene_tpm_mean = np.array(tpm.X.mean(axis=0)).reshape(-1)
    gene_tpm_stddev = var_sparse_matrix(tpm.X)**.5
    input_tpm_stats = pd.DataFrame([gene_tpm_mean, gene_tpm_stddev], index = ['_mean', '_std']).T
    save_df_to_npz(input_tpm_stats, pjoin(output_path, 'cNMF_split1/cnmf_tmp/cNMF_split1.tpm_stats.df.npz'))

    print("Computing highly variable genes and normalized counts on split1...")
    highvargenes = None
    genels_directory = pjoin(output_path, 'cNMF_split1/cNMF_split1.overdispersed_genes.txt')
    norm_counts = get_norm_counts(adata_split1, tpm, genels_directory, 
        num_highvar_genes=numgenes, high_variance_genes_filter=highvargenes)

    sc.write(pjoin(output_path, 'cNMF_split1/cnmf_tmp/cNMF_split1.norm_counts.h5ad'),norm_counts)

def check_dir_exists(path):
    """
    Checks if directory already exists or not and creates it if it doesn't
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def save_df_to_text(obj, filename):
    obj.to_csv(filename, sep='\t')
def save_df_to_npz(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)

def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj

def load_df_from_txt(filename):
    with open(filename) as f:
        obj = [int(line.rstrip('\n')) for line in f]
    return obj
def compute_tpm(input_counts):
    """
    Default TPM normalization
    input_counts: one split of raw counts
    """
    tpm = input_counts.copy()
    sc.pp.normalize_per_cell(tpm, counts_per_cell_after=1e6)
    return(tpm)

def var_sparse_matrix(X):
    mean = np.array(X.mean(axis=0)).reshape(-1)
    Xcopy = X.copy()
    Xcopy.data **= 2
    var = np.array(Xcopy.mean(axis=0)).reshape(-1) - (mean**2)
    return(var)

def get_highvar_genes_sparse(expression,  numgenes=None):
    # Find high variance genes within those cells
    gene_mean = np.array(expression.mean(axis=0)).astype(float).reshape(-1)
    E2 = expression.copy(); E2.data **= 2; gene2_mean = np.array(E2.mean(axis=0)).reshape(-1)
    gene_var = pd.Series(gene2_mean - (gene_mean**2))
    del(E2)
    gene_mean = pd.Series(gene_mean)
    gene_fano = gene_var / gene_mean
    # Find parameters for expected fano line
    top_genes = gene_mean.sort_values(ascending=False)[:20].index
    A = (np.sqrt(gene_var)/gene_mean)[top_genes].min()
    
    w_mean_low, w_mean_high = gene_mean.quantile([0.10, 0.90])
    w_fano_low, w_fano_high = gene_fano.quantile([0.10, 0.90])
    winsor_box = ((gene_fano > w_fano_low) &
                    (gene_fano < w_fano_high) &
                    (gene_mean > w_mean_low) &
                    (gene_mean < w_mean_high))
    fano_median = gene_fano[winsor_box].median()
    B = np.sqrt(fano_median)
    gene_expected_fano = (A**2)*gene_mean + (B**2)
    fano_ratio = (gene_fano/gene_expected_fano)

    # Identify high var genes
    if numgenes is not None:
        highvargenes = fano_ratio.sort_values(ascending=False).index[:numgenes]
        high_var_genes_ind = fano_ratio.index.isin(highvargenes)
        T=None

    else:
        if not expected_fano_threshold:
            #230227
            #T = (1. + gene_counts_fano[winsor_box].std())
            T = (1. + gene_fano[winsor_box].std())
        else:
            T = expected_fano_threshold

        #230227
        #high_var_genes_ind = (fano_ratio > T) & (gene_counts_mean > minimal_mean)
        high_var_genes_ind = (fano_ratio > T) & (gene_mean > minimal_mean)

    gene_counts_stats = pd.DataFrame({
        'mean': gene_mean,
        'var': gene_var,
        'high_var': high_var_genes_ind,
        })
    return(gene_counts_stats)

def get_norm_counts(input_counts, tpm, write_genels_directory, high_variance_genes_filter = None, num_highvar_genes = None):
    if high_variance_genes_filter is None:
        ## Get list of high-var genes if one wasn't provided
        if sparse.issparse(tpm.X):
            gene_counts_stats = get_highvar_genes_sparse(tpm.X, numgenes=num_highvar_genes)  
        else:
            raise TypeError("Input data has to be sparse!")
        high_variance_genes_filter = list(tpm.var.index[gene_counts_stats.high_var.values])
            
    ## Subset out high-variance genes
    norm_counts = input_counts[tpm.obs_names, high_variance_genes_filter]

    ## Scale genes to unit variance
    if sparse.issparse(tpm.X):
        sc.pp.scale(norm_counts, zero_center=False)
        if np.isnan(norm_counts.X.data).sum() > 0:
            print('Warning NaNs in normalized counts matrix')                       
    else:
        norm_counts.X /= norm_counts.X.std(axis=0, ddof=1)
        if np.isnan(norm_counts.X).sum().sum() > 0:
            print('Warning NaNs in normalized counts matrix')                    
    
    ## Save a \n-delimited list of the high-variance genes used for factorization
    open(write_genels_directory, 'w').write('\n'.join(high_variance_genes_filter))

    ## Check for any cells that have 0 counts of the overdispersed genes
    zerocells = norm_counts.X.sum(axis=1)==0
    if zerocells.sum()>0:
        examples = norm_counts.obs[zerocells].index
        print('Warning: %d cells have zero counts of overdispersed genes. E.g. %s' % (zerocells.sum(), examples[0]))
        print('Consensus step may not run when this is the case')
    
    return(norm_counts)

def fast_euclidean(mat):
    D = mat.dot(mat.T)
    squared_norms = np.diag(D).copy()
    D *= -2.0
    D += squared_norms.reshape((-1,1))
    D += squared_norms.reshape((1,-1))
    D[D<0]=0 # 030123 modified by XL
    D = np.sqrt(D)
    #D[D < 0] = 0
    return D
    
def fast_ols_all_cols(X, Y):
    pinv = np.linalg.pinv(X)
    beta = np.dot(pinv, Y)
    return(beta)

def fast_ols_all_cols_df(X,Y):
    beta = fast_ols_all_cols(X, Y)
    beta = pd.DataFrame(beta, index=X.columns, columns=Y.columns)
    return(beta)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= 'Input for running model')
    parser.add_argument('--output-dir', type=str, help='[all] Output directory. All output will be placed in [output-dir]/[name]/...', nargs='?', default='.')
    parser.add_argument('-c', '--counts', type=str, help='[prepare] Input (cell x gene) counts matrix as h5ad file')
    parser.add_argument('--numgenes', type=int, help='[prepare] Number of high variance genes to use for matrix factorization.', default=2000)
    args = parser.parse_args()


    raw_path = args.counts
    out_folder_path = args.output_dir
    counts_split(raw_path, out_folder_path, numgenes=args.numgenes)

    