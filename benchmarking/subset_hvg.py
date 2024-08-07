import numpy as np
import pandas as pd
import os, errno
import scipy.sparse as sp

from scipy.spatial.distance import squareform
from sklearn.decomposition import non_negative_factorization
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.utils import sparsefuncs



import matplotlib.pyplot as plt

import scanpy as sc




def save_df_to_npz(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)

def save_df_to_text(obj, filename):
    obj.to_csv(filename, sep='\t')

def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj




def get_highvar_genes_sparse(expression, expected_fano_threshold=None,
                       minimal_mean=0.5, numgenes=None):
    # Find high variance genes within those cells
    gene_mean = np.array(expression.mean(axis=0)).astype(float).reshape(-1)
    E2 = expression.copy(); E2 **= 2; gene2_mean = np.array(E2.mean(axis=0)).reshape(-1)
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

            T = (1. + gene_fano[winsor_box].std())
        else:
            T = expected_fano_threshold
      
        high_var_genes_ind = (fano_ratio > T) & (gene_mean > minimal_mean)

    gene_counts_stats = pd.DataFrame({
        'mean': gene_mean,
        'var': gene_var,
        'fano': gene_fano,
        'expected_fano': gene_expected_fano,
        'high_var': high_var_genes_ind,
        'fano_ratio': fano_ratio
        })
    gene_fano_parameters = {
            'A': A, 'B': B, 'T':T, 'minimal_mean': minimal_mean,
        }
    return(gene_counts_stats, gene_fano_parameters)




def compute_tpm(input_counts):
    """
    Default TPM normalization
    """
    tpm = input_counts.copy()
    sc.pp.normalize_per_cell(tpm, counts_per_cell_after=1e6, copy = True)
    return(tpm)


import sys, argparse
parser = argparse.ArgumentParser()

parser.add_argument('--counts', type=str, help='directory to raw count npz file', nargs='?', default='cNMF')
parser.add_argument('--num_highvar_genes', type=int, help='number of highly variable genes', nargs='?', default='.')

args = parser.parse_args()


input_counts = load_df_from_npz(args.counts)
# input_counts = t(input_counts) # Only GOSH needs this transpose transformation!!

tpm = compute_tpm(input_counts)
(gene_counts_stats, gene_fano_params) = get_highvar_genes_sparse(tpm, numgenes=args.num_highvar_genes)  

high_variance_genes_filter = list(tpm.columns.values[gene_counts_stats.high_var.values])
filtered_counts = input_counts[high_variance_genes_filter]

save_df_to_npz(filtered_counts, "benchmarking/filtered_counts_%d.npz" % (args.num_highvar_genes))
