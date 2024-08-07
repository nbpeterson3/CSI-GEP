#!/usr/bin/env python
# coding: utf-8

import os, sys

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse


def save_df_to_npz(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)


# hyper-parameters
ngenes=25000
ncells=250000
ngroups=40
libloc=7.64
libscale=0.78
mean_rate=7.68
mean_shape=0.34
expoutprob=0.00286
expoutloc=6.15
expoutscale=0.49
diffexpprob=0.025
diffexpdownprob=0.
diffexploc=1.0
diffexpscale=1.0
bcv_dispersion=0.448
bcv_dof=22.087
ndoublets=0
nproggenes=1000
progdownprob=0.
progdeloc=1.0
progdescale=1.0
progcellfrac=.35
proggroups=[1,2,3,4]

minprogusage=.1
maxprogusage=.7
seed=9458 
init_ncells = ncells + ndoublets
groupprob = [1/float(ngroups)]*ngroups

np.random.seed(seed)


# Get cell params

# Simulate groups

groupid = np.random.choice(np.arange(1, ngroups+1), size=init_ncells, p=groupprob)
groups = np.unique(groupid)
libsize = np.random.lognormal(mean=libloc, sigma=libscale, size=init_ncells)
cellnames = ['Cell%d' % i for i in range(1, init_ncells+1)]
cellparams = pd.DataFrame([groupid, libsize], index=['group', 'libsize'], columns=cellnames).T
cellparams['group'] = cellparams['group'].astype(int)

# Get gene params
basegenemean = np.random.gamma(shape=mean_shape,scale=1./mean_rate,size=ngenes)
is_outlier = np.random.choice([True, False], size=ngenes,
                                      p=[expoutprob,1-expoutprob])
outliers = np.random.lognormal(mean=expoutloc,sigma=expoutscale,size=is_outlier.sum())
outlier_ratio = np.ones(shape=ngenes)
outlier_ratio[is_outlier] = outliers
gene_mean = basegenemean.copy()
median = np.median(basegenemean)
gene_mean[is_outlier] = outliers*median
genenames = ['Gene%d' % i for i in range(1, ngenes+1)]
geneparams = pd.DataFrame([basegenemean, is_outlier, outlier_ratio, gene_mean], 
                          index=['BaseGeneMean', 'is_outlier', 'outlier_ratio', 'gene_mean'], columns=genenames).T
geneparams.iloc[is_outlier]


# Simulate programs
geneparams['prog_gene'] = False
proggenes = geneparams.index[-nproggenes:]
geneparams.loc[proggenes, 'prog_gene'] = True
DEratio = np.random.lognormal(mean=progdeloc,sigma=progdescale,size=nproggenes)
DEratio[DEratio<1] = 1 / DEratio[DEratio<1]
is_downregulated = np.random.choice([True, False],size=len(DEratio),
                                            p=[progdownprob,1-progdownprob])
DEratio[is_downregulated] = 1. / DEratio[is_downregulated]
all_DE_ratio = np.ones(ngenes)
all_DE_ratio[-nproggenes:] = DEratio[-nproggenes:]
prog_mean = geneparams['gene_mean']*all_DE_ratio
geneparams['prog_genemean'] = prog_mean

cellparams['has_program'] = False
cellparams.loc[:, 'program_usage'] = 0


for g in proggroups:
    groupcells = cellparams.index[cellparams['group']==g]
    hasprog = np.random.choice([True, False], size=len(groupcells),
                              p=[progcellfrac,1-progcellfrac])
    cellparams.loc[groupcells[hasprog], 'has_program'] = True
    usages = np.random.uniform(low=minprogusage,
                               high=maxprogusage,
                               size=len(groupcells[hasprog]))
    cellparams.loc[groupcells[hasprog], 'program_usage'] = usages
    


# Simulate group DE
proggene = geneparams['prog_gene'].values

for group in groups:
    isDE = np.random.choice([True, False], size=ngenes,
                              p=[diffexpprob,1-diffexpprob])
    isDE[proggene] = False # Program genes shouldn't be differentially expressed between groups
   
    DEratio = np.random.lognormal(mean=diffexploc,
                                  sigma=diffexpscale,
                                  size=isDE.sum())
    DEratio[DEratio<1] = 1 / DEratio[DEratio<1]
    is_downregulated = np.random.choice([True, False],
                                    size=len(DEratio),
                                    p=[diffexpdownprob,1-diffexpdownprob])
    DEratio[is_downregulated] = 1. / DEratio[is_downregulated]
    all_DE_ratio = np.ones(ngenes)
    all_DE_ratio[isDE] = DEratio
    group_mean = geneparams['gene_mean']*all_DE_ratio

    deratiocol = 'group%d_DEratio' % group
    groupmeancol = 'group%d_genemean' % group
    geneparams[deratiocol] = all_DE_ratio
    geneparams[groupmeancol] = group_mean


# Get cell gene means

group_genemean = geneparams.loc[:,[x for x in geneparams.columns if ('_genemean' in x) and ('group' in x)]].T.astype(float)
group_genemean = group_genemean.div(group_genemean.sum(axis=1), axis=0)
ind = cellparams['group'].apply(lambda x: 'group%d_genemean' % x)
noprogcells = cellparams['has_program']==False
hasprogcells = cellparams['has_program']==True
progcellmean = group_genemean.loc[ind[hasprogcells], :]
progcellmean.index = ind.index[hasprogcells]
progcellmean = progcellmean.multiply(1-cellparams.loc[hasprogcells, 'program_usage'], axis=0)


progmean = geneparams.loc[:,['prog_genemean']]
progmean = progmean.div(progmean.sum(axis=0), axis=1)
progusage = cellparams.loc[progcellmean.index, ['program_usage']]
progusage.columns = ['prog_genemean']
progcellmean += progusage.dot(progmean.T)
progcellmean = progcellmean.astype(float)


noprogcellmean = group_genemean.loc[ind[noprogcells],:]
noprogcellmean.index = ind.index[noprogcells]

cellgenemean = pd.concat([noprogcellmean, progcellmean], axis=0)

del(progcellmean, noprogcellmean)
cellgenemean = cellgenemean.reindex(index=cellparams.index)


normfac = (cellparams['libsize'] / cellgenemean.sum(axis=1)).values
for col in cellgenemean.columns:
    cellgenemean[col] = cellgenemean[col].values*normfac


# Adjust means bcv (biological coefficients of variation)

bcv = bcv_dispersion + (1 / np.sqrt(cellgenemean))
chisamp = np.random.chisquare(bcv_dof, size=ngenes)
bcv = bcv*np.sqrt(bcv_dof / chisamp)
updatedmean = np.random.gamma(shape=1/(bcv**2), scale=cellgenemean*(bcv**2))


bcv = pd.DataFrame(bcv, index=cellnames, columns=genenames)
updatedmean = pd.DataFrame(updatedmean, index=cellnames, columns=genenames)


# ## Simulate counts

counts = pd.DataFrame(np.random.poisson(lam=updatedmean),
                      index=cellnames, columns=genenames)



outdir = "simulated_data"
os.mkdir(outdir)
save_df_to_npz(cellparams, '%s/cellparams' % outdir)
save_df_to_npz(geneparams, '%s/geneparams' % outdir)
save_df_to_npz(counts, '%s/counts' % outdir)


adata = ad.AnnData(counts, dtype = 'float64')
sparse_X = sparse.csr_matrix(adata.X)
adata.X = sparse_X
adata.write_h5ad('%s/counts.h5ad' % outdir)

