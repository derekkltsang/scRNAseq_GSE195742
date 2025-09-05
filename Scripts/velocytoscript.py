#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 08:57:05 2021

@author: derek
"""
import matplotlib as plt
import numpy as np
import pandas as pd
import scvelo as scv
import anndata
import velocyto as vcy
import os
os.getcwd()
os.chdir('/home/derek/Documents/velocyto')


vlm = vcy.VelocytoLoom("S02.loom")

vlm.normalize("S", size=True, log=True)
vlm.S_norm  # contains log normalized

vlm.plot_fractions()

# Save
vlm.dump_hdf5("S02_velocyto_analysis")
# Load
load_velocyto_hdf5("my_velocyto_analysis.hdf5")

# ------------
import anndata
import scvelo as scv
import pandas as pd
from IPython.display import display
import numpy as np
import matplotlib as plt
%load_ext rpy2.ipython
# pip install rpy2

adata = anndata.read_loom("S03.loom")

sample_obs = pd.read_csv("SPF_mito_SEC_cellID_obs_v2.csv")
umap = pd.read_csv("SPF_mito_SEC_cell_embeddings_v2.csv")
clusters = pd.read_csv("SPF_mito_SEC_clusters_v2.csv")

#filter based on seurat IDs
adata = adata[np.isin(adata.obs.index,sample_obs["x"])]

#add UMAP coordinates
umap = pd.read_csv("SPF_mito_SEC_cell_embeddings_v2.csv")

adata.obs.index

adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {'CellID':'Cell ID'})
umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})

umap_ordered = adata_index.merge(umap, on = "Cell ID")

umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values

#add seurat clusters 
clusters = pd.read_csv("SPF_mito_SEC_clusters_v2.csv")
adata.obs.index

adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {'CellID':'Cell ID'})
clusters = clusters.rename(columns = {'Unnamed: 0':'Cell ID'})

clusters_ordered = adata_index.merge(clusters, on = "Cell ID")

clusters_ordered = clusters_ordered.iloc[:,1:]
adata.obs['clusters'] = clusters_ordered.values

#add seurat clusters colours
#ident_colours = ['#64B200','#00BD5C', '#00C1A7','#EF67EB','#AEA200','#00A6FF','#B385FF','#F8766D','#00BADE','#DB8E00','#FF63B6']
ident_colours = ['#EF67EB','#00A6FF', '#B385FF','#F8766D', '#00BADE', '#DB8E00', '#AEA200', '#FF63B6', '#64B200','#00BD5C','#00C1A7']

adata.uns['clusters_colors'] = ident_colours

#scVelo
adata.var_names_make_unique()
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
#scv.tl.recover_dynamics(adata, n_jobs=1) #only run if using dynamic model
#stochastic, deterministic, dynamical
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis = 'umap')


scv.pl.velocity_embedding_stream(adata, basis='umap', linewidth=1, legend_fontsize=8 ,legend_loc='lower right', size=400, figsize=(6,5), dpi=150,save=('SPF_scvelo_SEC_stream_stc_v2.png'))
scv.pl.proportions(adata)

adata.write_loom("S02_mito_SEC_scvelo",write_obsm_varm=True)

#check velocity of genes
scv.pl.velocity(adata, ['Hes1',  'Atoh1', 'Gfi1', 'Mki67', 'Alpi'], ncols=2,size=200, linewidth=2, fontsize=25, figsize=(10,8), dpi=150,save=('SPF_scvelo_stc_KeyGenes.png'))

#velocity genes
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv (r'/home/derek/Documents/velocyto/SPF_scVelo_SEC_stc_VeloGenes.csv', index = False, header=True)

kwargs = dict(frameon=False, fontsize=25, size=100, linewidth=3, figsize=(5,4), dpi=150)
#scv.pl.scatter(adata, df['0-SC/TA'][:5], ylabel='0-SC/TA', **kwargs)
scv.pl.scatter(adata, df['ISC'][:5], ylabel='ISC',**kwargs, save=('SPF_scvelo_SEC_stc_ISC.png'))
scv.pl.scatter(adata, df['TA1'][:5], ylabel='TA1',**kwargs, save=('SPF_scvelo_SEC_stc_TA1.png'))
scv.pl.scatter(adata, df['TA2'][:5], ylabel='TA2',**kwargs, save=('SPF_scvelo_SEC_stc_TA2.png'))
scv.pl.scatter(adata, df['V1-V2'][:5], ylabel='V1-V2',**kwargs, save=('SPF_scvelo_SEC_stc_V1-V2.png'))
scv.pl.scatter(adata, df['V3-V4'][:5], ylabel='V3-V4',**kwargs, save=('SPF_scvelo_SEC_stc_V3-V4.png'))
scv.pl.scatter(adata, df['V5-V6'][:5], ylabel='V5-V6',**kwargs, save=('SPF_scvelo_SEC_stc_V5-V6.png'))
scv.pl.scatter(adata, df['PC'][:5], ylabel='PC',**kwargs, save=('SPF_scvelo_stc_SEC_PC.png'))
scv.pl.scatter(adata, df['GC1'][:5], ylabel='GC1',**kwargs, save=('SPF_scvelo_SEC_stc_GC1.png'))
scv.pl.velocity(adata, df['GC2'][:5], ylabel='GC2',**kwargs, save=('SPF_scvelo_SEC_stc_GC2.png'))
scv.pl.scatter(adata, df['EE'][:5], ylabel='EE',**kwargs, save=('SPF_scvelo_stc_SEC_EE.png'))
scv.pl.scatter(adata, df['Tuft'][:5], ylabel='Tuft',**kwargs, save=('SPF_scvelo_SEC_stc_Tuft.png'))

#Speed & Coherence
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], size=100,figsize=(10,8), dpi=150, save=('SPF_scvelo_stc_speed.png'))

df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)
display(df)
df.to_csv (r'/home/derek/Documents/velocyto/SPF_scVelo_SEC_stc_Speed.csv', index = False, header=True)

scv.pl.velocity_graph(adata, threshold=.1, figsize=(10,8), size=150, dpi=150, save=('SPF_scvelo_SEC_stc_connections.png'))

#PAGA path
!pip install python-igraph --upgrade --quiet

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.velocity_confidence(adata)
scv.tl.paga(adata, groups='clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5)

