import pandas as pd
from os.path import join as pjoin
import os, errno, gc
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'
import argparse
import numpy as np
import scanpy as sc
from scipy import io, sparse
from random import choices

import torch
from torchnmf.nmf import NMF
from torchnmf.metrics import kl_div
print("Is CUDA available: ", torch.cuda.is_available())

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import squareform
from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list


def consensus(k, n_iter, output_tmp_path, output_path, density_threshold_str='0.5', local_neighborhood_size = 0.30):
    tpm = sc.read(output_tmp_path+'.tpm.h5ad')
    tpm_stats = load_df_from_npz(output_tmp_path+'.tpm_stats.df.npz')
    norm_counts = sc.read(output_tmp_path+'.norm_counts.h5ad')


    density_threshold_repl = density_threshold_str.replace('.', '_')
    density_threshold = float(density_threshold_str)
    n_neighbors = int(local_neighborhood_size * n_iter)

    combined_spectra = None
    # Combine all W matrices
    spectra_labels = []
    for i in range(n_iter):
        spectra = load_df_from_npz(output_tmp_path +'.spectra.k_%d.iter_%d.df.npz' % (k, i))
        if combined_spectra is None:
            combined_spectra = np.zeros((n_iter, k, spectra.shape[1]))
        combined_spectra[i, :, :] = spectra.values

        for t in range(k):
            spectra_labels.append('iter%d_topic%d'%(i, t+1))

    combined_spectra = combined_spectra.reshape(-1, combined_spectra.shape[-1])
    combined_spectra = pd.DataFrame(combined_spectra, columns=spectra.columns, index=spectra_labels)

    save_df_to_npz(combined_spectra, output_tmp_path +'.spectra.k_%d.merged.df.npz'%k)

    # Rescale topics such to length of 1.
    l2_spectra = (combined_spectra.T/np.sqrt((combined_spectra**2).sum(axis=1))).T

    topics_dist = None
    if os.path.isfile(output_tmp_path +'.local_density_cache.k_%d.merged.df.npz' % k):
        local_density = load_df_from_npz(output_tmp_path +'.local_density_cache.k_%d.merged.df.npz' % k)
    else:
        # first find the full distance matrix
        topics_dist = fast_euclidean(l2_spectra.values)
        #   partition based on the first n neighbors
        partitioning_order  = np.argpartition(topics_dist, n_neighbors+1)[:, :n_neighbors+1]
        #   find the mean over those n_neighbors (excluding self, which has a distance of 0)
        distance_to_nearest_neighbors = topics_dist[np.arange(topics_dist.shape[0])[:, None], partitioning_order]
        local_density = pd.DataFrame(distance_to_nearest_neighbors.sum(1)/(n_neighbors),
                                     columns=['local_density'],
                                     index=l2_spectra.index)
        save_df_to_npz(local_density, 
            output_tmp_path +'.local_density_cache.k_%d.merged.df.npz' % k)
        del(partitioning_order)
        del(distance_to_nearest_neighbors)

    density_filter = local_density.iloc[:, 0] < density_threshold
    l2_spectra = l2_spectra.loc[density_filter, :]

    kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1, algorithm = "elkan")
    kmeans_model.fit(l2_spectra)
    kmeans_cluster_labels = pd.Series(kmeans_model.labels_+1, index=l2_spectra.index)

    # Find median usage for each gene across cluster
    median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()

    # Normalize median spectra to probability distributions.
    median_spectra = (median_spectra.T/median_spectra.sum(1)).T

    # Compute the silhouette score
    stability = silhouette_score(l2_spectra.values, kmeans_cluster_labels, metric='euclidean')

    norm_counts.X = norm_counts.X.astype(np.float64)

    nnz_i, nnz_j, v = sparse.find(norm_counts.X)
    i = torch.from_numpy(np.array([nnz_i, nnz_j]))
    if torch.cuda.is_available():
        S = torch.sparse_coo_tensor(i, v,
                                 dtype=torch.float,
                                 device=torch.device('cuda:0'))
        m = NMF(rank = k, W = torch.tensor(median_spectra.transpose().values), trainable_W = False, trainable_H = True, H =  torch.rand(norm_counts.shape[0], k)).cuda()
        # beta = 1 - KL divergence; beta = 2 - Euclidean
        m.fit(S,verbose=False, beta=2, max_iter=1000, tol = 1e-6)
        usages = pd.DataFrame(list(m.H.detach().cpu().numpy()), norm_counts.obs.index, columns=median_spectra.index)
    else:
        print("GPU not available, calculating consensus matrix on CPU...")
        S = torch.sparse_coo_tensor(i, v,
                                 dtype=torch.float)
        m = NMF(rank = k, W = torch.tensor(median_spectra.transpose().values), trainable_W = False, trainable_H = True, H =  torch.rand(norm_counts.shape[0], k))
        # beta = 1 - KL divergence; beta = 2 - Euclidean
        m.fit(S,verbose=False, beta=2, max_iter=1000, tol = 1e-6)
        usages = pd.DataFrame(list(m.H.detach().numpy()), norm_counts.obs.index, columns=median_spectra.index)
 
    rf_pred_norm_counts = usages.dot(median_spectra)

    # Compute prediction error as a frobenius norm
    if sparse.issparse(norm_counts.X):
        prediction_error = ((norm_counts.X.todense() - rf_pred_norm_counts)**2).sum().sum()
    else:
        prediction_error = ((norm_counts.X - rf_pred_norm_counts)**2).sum().sum()

    consensus_stats = pd.DataFrame([k, density_threshold, stability, prediction_error],
                index = ['k', 'local_density_threshold', 'stability', 'prediction_error'],
                columns = ['stats'])

    print(consensus_stats)

    save_df_to_npz(median_spectra, output_tmp_path + '.spectra.k_%d.dt_%s.consensus.df.npz'%(k, density_threshold_repl))
    save_df_to_npz(usages, output_tmp_path + '.usages.k_%d.dt_%s.consensus.df.npz'%(k, density_threshold_repl))
    save_df_to_npz(consensus_stats, output_tmp_path + '.stats.k_%d.dt_%s.df.npz'%(k, density_threshold_repl))
    save_df_to_text(median_spectra, output_path + '.spectra.k_%d.dt_%s.consensus.txt'%(k, density_threshold_repl))
    save_df_to_text(usages, output_path + '.usages.k_%d.dt_%s.consensus.txt'%(k, density_threshold_repl))
    
    if sparse.issparse(tpm.X):
        norm_tpm = (np.array(tpm.X.todense()) - tpm_stats['_mean'].values) / (tpm_stats['_std'].values +1e-10)
    else:
        norm_tpm = (tpm.X - tpm_stats['_mean'].values) / (tpm_stats['_std'].values + 1e-10)
    
    if norm_tpm.dtype != np.float64:
        norm_tpm = norm_tpm.astype(np.float64)
    
    usage_coef = fast_ols_all_cols(usages.values, norm_tpm)
    usage_coef = pd.DataFrame(usage_coef, index=usages.columns, columns=tpm.var.index)

    save_df_to_npz(usage_coef, output_path + '.gene_spectra_score.k_%d.dt_%s.df.npz'%(k, density_threshold_repl))
    save_df_to_text(usage_coef, output_path + '.gene_spectra_score.k_%d.dt_%s.txt'%(k, density_threshold_repl))


def factorize(count_path, output_path, k, n_iter, verbose = False):
    if torch.cuda.is_available():
        norm_counts = sc.read(count_path)
        nnz_i, nnz_j, v = sparse.find(norm_counts.X)
        i = torch.from_numpy(np.array([nnz_i, nnz_j]))
        S = torch.sparse_coo_tensor(i, v,
                             dtype=torch.float,
                             device=torch.device('cuda:0'))
        del nnz_i, nnz_j, i, v
        torch.cuda.empty_cache()

        for i in range(n_iter):
            torch.cuda.empty_cache()
            m = NMF(S.shape, rank = k).cuda()
            # beta = 1 - KL divergence; beta = 2 - Euclidean
            m.fit(S,verbose=verbose, beta=2, max_iter=500, tol = 1e-6)
            spectra = pd.DataFrame(list(m.W.detach().cpu().numpy()), index= norm_counts.var_names, columns = None).transpose()
            np.savez_compressed(output_path+'.spectra.k_%d.iter_%d.df.npz' % (k, i), 
                data = spectra.values, index = spectra.index.values, columns = spectra.columns.values)
            torch.cuda.empty_cache()
            gc.collect()
            del m, spectra

    else:
        print("GPU not found, factorizing on CPU...")
        norm_counts = sc.read(count_path)
        nnz_i, nnz_j, v = sparse.find(norm_counts.X)
        i = torch.from_numpy(np.array([nnz_i, nnz_j]))
        S = torch.sparse_coo_tensor(i, v,
                             dtype=torch.float)
        del nnz_i, nnz_j, i, v

        for i in range(n_iter):
            m = NMF(S.shape, rank = k)
            # beta = 1 - KL divergence; beta = 2 - Euclidean
            m.fit(S,verbose = verbose, beta=2, max_iter=500, tol = 1e-6)
            spectra = pd.DataFrame(list(m.W.detach().numpy()), index= norm_counts.var_names, columns = None).transpose()
            np.savez_compressed(output_path+'.spectra.k_%d.iter_%d.df.npz' % (k, i), 
                data = spectra.values, index = spectra.index.values, columns = spectra.columns.values)
            gc.collect()
            del m, spectra




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
    parser.add_argument('-k', '--components', type=int, help='[factorize] Numper of components (k) for matrix factorization. Input as list')
    parser.add_argument('-n', '--n-iter', type=int, help='[factorize] Numper of factorization replicates', default=100)
    args = parser.parse_args()


    out_folder_path = args.output_dir
    k = args.components
    niter = args.n_iter
    
    # split0 
    count_path = pjoin(out_folder_path, 'cNMF_split0/cnmf_tmp/cNMF_split0.norm_counts.h5ad') # directory of split0 normalized data
    out_tmp_path = pjoin(out_folder_path, 'cNMF_split0/cnmf_tmp/cNMF_split0') # directory to where all tmp files (from factorize and combine) write
    out_path = pjoin(out_folder_path, 'cNMF_split0/cNMF_split0') # directory where full usages and spectra matrices store at


    factorize(count_path, out_tmp_path, k, niter, verbose = False)
    density_thres_str='2.0'
    consensus(k, niter, out_tmp_path, out_path, density_threshold_str=density_thres_str, local_neighborhood_size = 0.30)
    density_thres_str='0.1'
    consensus(k, niter, out_tmp_path, out_path, density_threshold_str=density_thres_str, local_neighborhood_size = 0.30)


    # split1 
    count_path = pjoin(out_folder_path, 'cNMF_split1/cnmf_tmp/cNMF_split1.norm_counts.h5ad') # directory of split0 normalized data
    out_tmp_path = pjoin(out_folder_path, 'cNMF_split1/cnmf_tmp/cNMF_split1') # directory to where all tmp files (from factorize and combine) write
    out_path = pjoin(out_folder_path, 'cNMF_split1/cNMF_split1') # directory where full usages and spectra matrices store at

    factorize(count_path, out_tmp_path, k, niter, verbose = False)
    density_thres_str='2.0'
    consensus(k, niter, out_tmp_path, out_path, density_threshold_str=density_thres_str, local_neighborhood_size = 0.30)
    density_thres_str='0.1'
    consensus(k, niter, out_tmp_path, out_path, density_threshold_str=density_thres_str, local_neighborhood_size = 0.30)
