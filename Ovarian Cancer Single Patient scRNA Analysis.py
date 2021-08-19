#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
import scanpy as sc
import numpy as np
import pandas as pd
import csv
from csv import reader
import itertools
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import pickle
import gseapy as gp
from gseapy.plot import gseaplot


# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[3]:


I = pickle.load(open("intron_network.pkl","rb"))
U = pickle.load(open("utr_network.pkl", "rb"))
N = pickle.load(open("non_network.pkl", "rb"))
S = pickle.load(open("similar_network.pkl", "rb"))


# # SPECTRUM-OV-009 Patient Analysis

# In[4]:


adata = sc.read_10x_mtx('SPECTRUM_OV_009/', cache=True)


# In[5]:


metadata = pd.read_csv("SPECTRUM-OV-009.tsv",sep=",")
celltypes = dict(zip(metadata["barcode"],metadata["cell_type"]))
ordered_celltypes = []
for barcode in adata.obs.index:
    ordered_celltypes.append(celltypes[barcode])
adata.obs["cell_type"] = ordered_celltypes
samples = dict(zip(metadata["barcode"],metadata["sample"]))
ordered_samples = []
for barcode in adata.obs.index:
    ordered_samples.append(samples[barcode])
adata.obs["sample"] = ordered_samples
adata.var_names_make_unique()


# In[6]:


adata.obs


# In[7]:


sc.pl.highest_expr_genes(adata, n_top=20, )


# In[8]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[9]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[10]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)


# In[11]:


sc.pp.calculate_qc_metrics(adata)


# In[12]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[13]:


adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)


# In[14]:


sc.pl.umap(adata,color="cell_type")
adata = adata[adata.obs["cell_type"]=="Ovarian.cancer.cell"]
sc.pl.umap(adata,color="cell_type")


# In[15]:


genes = adata.var_names
gene_list = genes.tolist()
adata.obs


# In[16]:


with open('RBP2GO_RBP_list.csv', newline='', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)
    RBP_list = list(itertools.chain.from_iterable(data))
valid_rbps = set(RBP_list).intersection(set(gene_list))


# In[17]:


sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)


# In[18]:


sc.tl.leiden(adata, resolution=0.2)
sc.pl.umap(adata, color=['leiden'])
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[19]:


sc.pl.rank_genes_groups_heatmap(adata, n_genes=50, use_raw=False, swap_axes=True, show_gene_labels=False,
                                vmin=-3, vmax=3, cmap='bwr')


# In[20]:


sc.pl.rank_genes_groups_tracksplot(adata, n_genes=5)


# In[21]:


valid_rbps = set()
for x in I.edges():
    valid_rbps.add(x[0])


# In[22]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for i, cluster in enumerate(degs):
    cluster_degs = degs[cluster].tolist()[:25]
    up_regulated_targets = set(cluster_degs).intersection(I.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(I.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    if not list(upstream_rbps) + list(up_regulated_targets) == []:
        enr = gp.enrichr(gene_list=list(upstream_rbps) + list(up_regulated_targets),
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        print(enr.results.head(50))
    subI = I.edge_subgraph(valid_edges)
    figure(figsize=(16, 8))
    ax = plt.subplot(1,2,1)
    nx.draw(subI, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["leiden"] == str(i)]
    sc.pl.umap(cluster_sub_umap,color="leiden",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[23]:


valid_rbps = set()
for x in U.edges():
    valid_rbps.add(x[0])


# In[24]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for i, cluster in enumerate(degs):
    cluster_degs = degs[cluster].tolist()[:25]
    up_regulated_targets = set(cluster_degs).intersection(U.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(U.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    if not list(upstream_rbps) + list(up_regulated_targets) == []:
        enr = gp.enrichr(gene_list=list(upstream_rbps) + list(up_regulated_targets),
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        print(enr.results.head(50))
    subU = U.edge_subgraph(valid_edges)
    figure(figsize=(16, 8))
    ax = plt.subplot(1,2,1)
    nx.draw(subU, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["leiden"] == str(i)]
    sc.pl.umap(cluster_sub_umap,color="leiden",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[25]:


valid_rbps = set()
for x in N.edges():
    valid_rbps.add(x[0])


# In[26]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for i, cluster in enumerate(degs):
    cluster_degs = degs[cluster].tolist()[:25]
    up_regulated_targets = set(cluster_degs).intersection(N.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(N.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    if not list(upstream_rbps) + list(up_regulated_targets) == []:
        enr = gp.enrichr(gene_list=list(upstream_rbps) + list(up_regulated_targets),
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        print(enr.results.head(50))
    subN = N.edge_subgraph(valid_edges)
    figure(figsize=(16, 8))
    ax = plt.subplot(1,2,1)
    nx.draw(subN, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["leiden"] == str(i)]
    sc.pl.umap(cluster_sub_umap,color="leiden",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[27]:


valid_rbps = set()
for x in S.edges():
    valid_rbps.add(x[0])


# In[28]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for i, cluster in enumerate(degs):
    cluster_degs = degs[cluster].tolist()[:25]
    up_regulated_targets = set(cluster_degs).intersection(S.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(S.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    if not list(upstream_rbps) + list(up_regulated_targets) == []:
        enr = gp.enrichr(gene_list=list(upstream_rbps) + list(up_regulated_targets),
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        print(enr.results.head(50))
    subS = S.edge_subgraph(valid_edges)
    figure(figsize=(16, 8))
    ax = plt.subplot(1,2,1)
    nx.draw(subS, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["leiden"] == str(i)]
    sc.pl.umap(cluster_sub_umap,color="leiden",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

