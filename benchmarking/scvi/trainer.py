import os
import torch
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import sys, argparse

def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj

parser = argparse.ArgumentParser()
parser.add_argument('--counts', type=str, help='directory to raw count npz file', nargs='?', default='cNMF')
parser.add_argument('--K0', type = int, help = 'dimensions of latent space', nargs = '?', default = 10)
parser.add_argument('--layers', type = int, help = 'number of fully connected layers', nargs = '?', default = 1)
parser.add_argument('--option', type = int, help = 'mode of fitting model', default = 1)
args = parser.parse_args()

input_counts = load_df_from_npz(args.counts)
adata = sc.AnnData(input_counts, 
    input_counts.index.to_frame(), 
    input_counts.columns.to_frame())

adata.layers["counts"] = adata.X.copy()  # preserve counts
mode = args.option

if mode == 1:
    ## search for K
    K = ' '.join([str(i) for i in range(2,61)] + [str(i) for i in range(70,210,10)])
    K_candidate  = [int(i) for i in K.split()]
    for x in K_candidate:
        scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts")
        model = scvi.model.SCVI(adata, gene_likelihood = "nb", latent_distribution = "ln", n_layers=args.layers, n_latent=x)
        model.train(max_epochs = 500, use_gpu = 'cuda:0', enable_progress_bar = False, early_stopping = True)
        print(f"k = {x}, marginal log-likelihood = {model.get_marginal_ll()}")
if mode == 2:
    ## fixed K
    K = args.K0
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts"
    )
    model = scvi.model.SCVI(adata, gene_likelihood = "nb", latent_distribution = "ln", n_layers=args.layers, n_latent=K)
    model.train(max_epochs = 500, use_gpu = 'cuda:0', enable_progress_bar = False)
    rep = model.get_latent_representation()
    np.savetxt("benchmarking/scvi/output/scvi_representation_%d.txt" % (K), rep)

