#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
import pickle
from scipy.sparse import csr_matrix
import scanpy as sc
import numpy as np
import pandas as pd
import csv
from csv import reader
import itertools
import glob
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import gseapy as gp
import pygraphviz
import collections


# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[3]:


adata = sc.read('spectrum_cancer_signatures.h5ad', cache=True)


# In[4]:


sc.pl.umap(adata,color="patient")
sc.pl.umap(adata,color="signature")


# In[5]:


sc.tl.rank_genes_groups(adata, 'signature', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# # Intron Cancer Patient-Network Analysis

# In[6]:


I = pickle.load(open("intron_network.pkl", "rb"))
valid_rbps = set()
for x in I.edges():
    valid_rbps.add(x[0])


# In[7]:


column = "signature"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:400]
    up_regulated_targets = set(cluster_degs).intersection(I.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(I.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subI = I.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subI, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subI, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[8]:


up_targets_all = list(set(up_targets_all))
sc.pl.matrixplot(adata, up_targets_all, 'patient',swap_axes=True, dendrogram=True, cmap='plasma', standard_scale="var")
sc.pl.matrixplot(adata, up_targets_all, 'signature',swap_axes=True, dendrogram=True, cmap='plasma', standard_scale="var")


# In[9]:


up_rbps_all = list(set(up_rbs_all).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, up_rbps_all, 'patient',swap_axes=True, dendrogram=True, cmap='Reds',standard_scale='var')
sc.pl.matrixplot(adata, up_rbps_all, 'signature',swap_axes=True, dendrogram=True, cmap='Reds',standard_scale='var')


# In[10]:


marker_dict = dict()
for r in up_rbps_all:
    xtargets = [x[1] for x in I.out_edges(r)]
    marker_dict[r] = list(set(xtargets).intersection(set(adata.var.index)))
sc.pl.matrixplot(adata, marker_dict, 'signature',swap_axes=False, dendrogram=True, cmap='cividis',standard_scale='var')


# In[11]:


adata.var.index.intersection(set(valid_rbps))


# In[12]:


pathway_markers = dict()
for pathway, pgenes in up_pathway_genes.items():
    pathway_markers = list(set(pgenes).intersection(set(adata.var.index)))
    if pathway_markers == []: continue
    sc.pl.matrixplot(adata, pathway_markers, 'signature', title=pathway,swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[13]:


for rbp in valid_rbps:
    targets = [x[1] for x in I.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, targets, 'patient', title=rbp,swap_axes=True, dendrogram=True)


# In[14]:


for rbp in set(valid_rbps).intersection(set(adata.var.index.tolist())):
    targets = [x[1] for x in I.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, [rbp] + targets, 'signature', title=rbp,swap_axes=True, dendrogram=True)


# In[15]:


vrbs = list(set(valid_rbps).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, vrbs, 'patient', swap_axes=True, dendrogram=True, cmap='Reds')
sc.pl.matrixplot(adata, vrbs, 'signature',swap_axes=True, dendrogram=True, cmap='Reds',standard_scale='var')


# In[16]:


zdata = adata[adata.obs["signature"]=="FBI"]
sc.pl.matrixplot(zdata, vrbs, 'patient', dendrogram=False, cmap='Blues')
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=False, cmap='Greens', standard_scale="var")


# In[17]:


xtargets = [x[1] for x in I.out_edges("RBM38")]
xtargets = ["RBM38"] + xtargets
z = adata.var.index.intersection(set(xtargets))
sc.pl.matrixplot(zdata, z, 'patient', dendrogram=False, cmap='Blues', standard_scale="var")


# # Untranslated Region Cancer Patient-Network Analysis

# In[18]:


U = pickle.load(open("utr_network.pkl","rb"))
valid_rbps = set()
for x in U.edges():
    valid_rbps.add(x[0])


# In[19]:


column = "signature"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:400]
    up_regulated_targets = set(cluster_degs).intersection(U.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(U.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subU = U.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subU, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subU, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[20]:


up_targets_all = list(set(up_targets_all))
sc.pl.matrixplot(adata, up_targets_all, 'patient',swap_axes=True, dendrogram=True, cmap='Blues', standard_scale="var")
sc.pl.matrixplot(adata, up_targets_all, 'signature',swap_axes=True, dendrogram=True, cmap='Blues', standard_scale="var")


# In[21]:


up_rbps_all = list(set(up_rbs_all).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, up_rbps_all, 'patient',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')
sc.pl.matrixplot(adata, up_rbps_all, 'signature',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[22]:


marker_dict = dict()
for r in up_rbps_all:
    xtargets = [x[1] for x in U.out_edges(r)]
    marker_dict[r] = list(set(xtargets).intersection(set(adata.var.index)))
sc.pl.matrixplot(adata, marker_dict, 'signature',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[23]:


adata.var.index.intersection(set(valid_rbps))


# In[24]:


pathway_markers = dict()
for pathway, pgenes in up_pathway_genes.items():
    pathway_markers = list(set(pgenes).intersection(set(adata.var.index)))
    if pathway_markers == []: continue
    sc.pl.matrixplot(adata, pathway_markers, 'signature', title=pathway,swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[25]:


for rbp in valid_rbps:
    targets = [x[1] for x in U.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, targets, 'patient', title=rbp,swap_axes=True, dendrogram=True)


# In[26]:


for rbp in set(valid_rbps).intersection(set(adata.var.index.tolist())):
    targets = [x[1] for x in U.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, [rbp] + targets, 'signature', title=rbp,swap_axes=True, dendrogram=True)


# In[27]:


vrbs = list(set(valid_rbps).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=True, cmap='Blues')


# In[28]:


zdata = adata[adata.obs["signature"]=="FBI"]
sc.pl.matrixplot(zdata, vrbs, 'patient', dendrogram=False, cmap='Blues')
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=False, cmap='Greens', standard_scale="var")


# In[29]:


xtargets = [x[1] for x in U.out_edges("KHDRBS3")]
xtargets = ['KHDRBS3'] + xtargets
z = adata.var.index.intersection(set(xtargets))
sc.pl.matrixplot(zdata, z, 'patient', dendrogram=False, cmap='Blues', standard_scale="var")


# # Non Coding RNA Cancer Patient-Network Analysis

# In[30]:


N = pickle.load(open("non_network.pkl", "rb"))
valid_rbps = set()
for x in N.edges():
    valid_rbps.add(x[0])


# In[31]:


column = "signature"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:400]
    up_regulated_targets = set(cluster_degs).intersection(N.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(N.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subN = N.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subN, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subN, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[32]:


up_targets_all = list(set(up_targets_all))
sc.pl.matrixplot(adata, up_targets_all, 'patient',swap_axes=True, dendrogram=True, cmap='Blues', standard_scale="var")
sc.pl.matrixplot(adata, up_targets_all, 'signature',swap_axes=True, dendrogram=True, cmap='Blues', standard_scale="var")


# In[33]:


up_rbps_all = list(set(up_rbs_all).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, up_rbps_all, 'patient',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')
sc.pl.matrixplot(adata, up_rbps_all, 'signature',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[34]:


marker_dict = dict()
for r in up_rbps_all:
    xtargets = [x[1] for x in N.out_edges(r)]
    marker_dict[r] = list(set(xtargets).intersection(set(adata.var.index)))
sc.pl.matrixplot(adata, marker_dict, 'signature',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[35]:


adata.var.index.intersection(set(valid_rbps))


# In[36]:


pathway_markers = dict()
for pathway, pgenes in up_pathway_genes.items():
    pathway_markers = list(set(pgenes).intersection(set(adata.var.index)))
    if pathway_markers == []: continue
    sc.pl.matrixplot(adata, pathway_markers, 'signature', title=pathway,swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[37]:


for rbp in valid_rbps:
    targets = [x[1] for x in N.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, targets, 'patient', title=rbp,swap_axes=True, dendrogram=True)


# In[38]:


for rbp in set(valid_rbps).intersection(set(adata.var.index.tolist())):
    targets = [x[1] for x in N.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, [rbp] + targets, 'signature', title=rbp,swap_axes=True, dendrogram=True)


# In[39]:


vrbs = list(set(valid_rbps).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=True, cmap='Blues')


# In[40]:


zdata = adata[adata.obs["signature"]=="FBI"]
sc.pl.matrixplot(zdata, vrbs, 'patient', dendrogram=False, cmap='Blues')
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=False, cmap='Greens', standard_scale="var")


# In[41]:


xtargets = [x[1] for x in N.out_edges("KHDRBS3")]
xtargets = ["KHDRBS3"] + xtargets
z = adata.var.index.intersection(set(xtargets))
sc.pl.matrixplot(zdata, z, 'patient', dendrogram=False, cmap='Blues', standard_scale="var")


# # Intersection Category Cancer Patient-Network Analysis

# In[42]:


S = pickle.load(open("similar_network.pkl", "rb"))
valid_rbps = set()
for x in S.edges():
    valid_rbps.add(x[0])


# In[43]:


column = "signature"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:400]
    up_regulated_targets = set(cluster_degs).intersection(S.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(S.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subS = S.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subS, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subS, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[44]:


up_targets_all = list(set(up_targets_all))
sc.pl.matrixplot(adata, up_targets_all, 'patient',swap_axes=True, dendrogram=True, cmap='Blues', standard_scale="var")
sc.pl.matrixplot(adata, up_targets_all, 'signature',swap_axes=True, dendrogram=True, cmap='Blues', standard_scale="var")


# In[45]:


up_rbps_all = list(set(up_rbs_all).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, up_rbps_all, 'patient',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')
sc.pl.matrixplot(adata, up_rbps_all, 'signature',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[46]:


marker_dict = dict()
for r in up_rbps_all:
    xtargets = [x[1] for x in S.out_edges(r)]
    marker_dict[r] = list(set(xtargets).intersection(set(adata.var.index)))
sc.pl.matrixplot(adata, marker_dict, 'signature',swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[47]:


adata.var.index.intersection(set(valid_rbps))


# In[48]:


pathway_markers = dict()
for pathway, pgenes in up_pathway_genes.items():
    pathway_markers = list(set(pgenes).intersection(set(adata.var.index)))
    if pathway_markers == []: continue
    sc.pl.matrixplot(adata, pathway_markers, 'signature', title=pathway,swap_axes=True, dendrogram=True, cmap='Blues',standard_scale='var')


# In[49]:


for rbp in valid_rbps:
    targets = [x[1] for x in S.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, targets, 'patient', title=rbp,swap_axes=True, dendrogram=True)


# In[50]:


for rbp in set(valid_rbps).intersection(set(adata.var.index.tolist())):
    targets = [x[1] for x in S.out_edges(rbp)]
    targets = list(set(targets).intersection(set(adata.var.index.tolist())))
    if len(targets) == 0: continue
    sc.pl.matrixplot(adata, [rbp] + targets, 'signature', title=rbp,swap_axes=True, dendrogram=True)


# In[51]:


vrbs = list(set(valid_rbps).intersection(set(adata.var.index.tolist())))
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=True, cmap='Blues')


# In[52]:


zdata = adata[adata.obs["signature"]=="FBI"]
sc.pl.matrixplot(zdata, vrbs, 'patient', dendrogram=False, cmap='Blues')
sc.pl.matrixplot(adata, vrbs, 'patient', dendrogram=False, cmap='Greens', standard_scale="var")


# In[53]:


xtargets = [x[1] for x in S.out_edges("RBM38")]
xtargets = ["RBM38"] + xtargets
z = adata.var.index.intersection(set(xtargets)) 
sc.pl.matrixplot(zdata, z, 'patient', dendrogram=False, cmap='Blues', standard_scale="var")


# # Control Patient-Network Analysis

# In[54]:


adata = sc.read('cancer_normal.h5ad', cache=True)


# In[55]:


sc.pl.umap(adata,color="cell_type")


# In[56]:


sc.tl.rank_genes_groups(adata, 'cell_type', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[57]:


I = pickle.load(open("intron_network.pkl", "rb"))
valid_rbps = set()
for x in I.edges():
    valid_rbps.add(x[0])


# In[58]:


column = "cell_type"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_targets_all_dict = dict()
up_rbps_all_dict = dict()
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:200]
    up_regulated_targets = set(cluster_degs).intersection(I.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(I.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all_dict[cluster] = list(up_regulated_targets)
    up_rbps_all_dict[cluster] = list(upstream_rbps)
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subI = I.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subI, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subI, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[59]:


cancer_specific_targets = set(up_targets_all_dict["Ovarian.cancer.cell"]).difference(set(up_targets_all_dict["Normal"]))
cancer_specific_rbps = set(up_rbps_all_dict["Ovarian.cancer.cell"]).difference(set(up_rbps_all_dict["Normal"]))
normal_specific_targets = set(up_targets_all_dict["Normal"]).difference(set(up_targets_all_dict["Ovarian.cancer.cell"]))
normal_specific_rbps = set(up_rbps_all_dict["Normal"]).difference(set(up_rbps_all_dict["Ovarian.cancer.cell"]))
enr = gp.enrichr(gene_list=list(cancer_specific_rbps) + list(upstream_rbps),
                     gene_sets=['MSigDB_Hallmark_2020'],
                     organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                     description='enrichr',
                     outdir='test/enrichr_kegg',
                     cutoff=0.05
                     )
print(enr.results.head(200))


# In[60]:


U = pickle.load(open("utr_network.pkl", "rb"))
valid_rbps = set()
for x in U.edges():
    valid_rbps.add(x[0])


# In[61]:


column = "cell_type"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_targets_all_dict = dict()
up_rbps_all_dict = dict()
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:200]
    up_regulated_targets = set(cluster_degs).intersection(U.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(U.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all_dict[cluster] = list(up_regulated_targets)
    up_rbps_all_dict[cluster] = list(upstream_rbps)
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subU = U.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subU, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subU, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[62]:


cancer_specific_targets = set(up_targets_all_dict["Ovarian.cancer.cell"]).difference(set(up_targets_all_dict["Normal"]))
cancer_specific_rbps = set(up_rbps_all_dict["Ovarian.cancer.cell"]).difference(set(up_rbps_all_dict["Normal"]))
normal_specific_targets = set(up_targets_all_dict["Normal"]).difference(set(up_targets_all_dict["Ovarian.cancer.cell"]))
normal_specific_rbps = set(up_rbps_all_dict["Normal"]).difference(set(up_rbps_all_dict["Ovarian.cancer.cell"]))
enr = gp.enrichr(gene_list=list(cancer_specific_rbps) + list(upstream_rbps),
                     gene_sets=['MSigDB_Hallmark_2020'],
                     organism='Human',
                     description='enrichr',
                     outdir='test/enrichr_kegg',
                     cutoff=0.05
                     )
print(enr.results.head(200))


# In[63]:


N = pickle.load(open("non_network.pkl", "rb"))
valid_rbps = set()
for x in N.edges():
    valid_rbps.add(x[0])


# In[64]:


column = "cell_type"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_targets_all_dict = dict()
up_rbps_all_dict = dict()
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:200]
    up_regulated_targets = set(cluster_degs).intersection(N.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(N.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all_dict[cluster] = list(up_regulated_targets)
    up_rbps_all_dict[cluster] = list(upstream_rbps)
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subN = N.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subN, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subN, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[65]:


cancer_specific_targets = set(up_targets_all_dict["Ovarian.cancer.cell"]).difference(set(up_targets_all_dict["Normal"]))
cancer_specific_rbps = set(up_rbps_all_dict["Ovarian.cancer.cell"]).difference(set(up_rbps_all_dict["Normal"]))
normal_specific_targets = set(up_targets_all_dict["Normal"]).difference(set(up_targets_all_dict["Ovarian.cancer.cell"]))
normal_specific_rbps = set(up_rbps_all_dict["Normal"]).difference(set(up_rbps_all_dict["Ovarian.cancer.cell"]))
enr = gp.enrichr(gene_list=list(cancer_specific_rbps) + list(upstream_rbps),
                     gene_sets=['MSigDB_Hallmark_2020'],
                     organism='Human',
                     description='enrichr',
                     outdir='test/enrichr_kegg',
                     cutoff=0.05
                     )
print(enr.results.head(200))


# In[66]:


S = pickle.load(open("similar_network.pkl", "rb"))
valid_rbps = set()
for x in S.edges():
    valid_rbps.add(x[0])


# In[67]:


column = "cell_type"
degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
up_targets_all = []
up_rbs_all = []
up_targets_all_dict = dict()
up_rbps_all_dict = dict()
up_pathway_genes = collections.defaultdict(list)
for cluster in degs:
    cluster_degs = degs[cluster].tolist()[:200]
    up_regulated_targets = set(cluster_degs).intersection(S.nodes())
    valid_edges = []
    upstream_rbps = set()
    for target in up_regulated_targets:
        new_edges = list(S.in_edges(target))
        valid_edges += new_edges
        for nedge in new_edges:
            upstream_rbps.add(nedge[0])
    print(cluster, up_regulated_targets, upstream_rbps)
    try:
        enr = gp.enrichr(gene_list=list(up_regulated_targets) + list(upstream_rbps),
                         gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='enrichr',
                         outdir='test/enrichr_kegg',
                         cutoff=0.05
                        )
        for n,g in zip(enr.results["Term"],enr.results["Genes"]):
            up_pathway_genes[n] += str(g).split(";")
    except Exception as e:
        continue
    up_targets_all_dict[cluster] = list(up_regulated_targets)
    up_rbps_all_dict[cluster] = list(upstream_rbps)
    up_targets_all += list(up_regulated_targets)
    up_rbs_all += list(upstream_rbps)
    subS = S.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subS, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subS, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs[column] == str(cluster)]
    sc.pl.umap(cluster_sub_umap,color=column,ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[68]:


cancer_specific_targets = set(up_targets_all_dict["Ovarian.cancer.cell"]).difference(set(up_targets_all_dict["Normal"]))
cancer_specific_rbps = set(up_rbps_all_dict["Ovarian.cancer.cell"]).difference(set(up_rbps_all_dict["Normal"]))
normal_specific_targets = set(up_targets_all_dict["Normal"]).difference(set(up_targets_all_dict["Ovarian.cancer.cell"]))
normal_specific_rbps = set(up_rbps_all_dict["Normal"]).difference(set(up_rbps_all_dict["Ovarian.cancer.cell"]))
enr = gp.enrichr(gene_list=list(cancer_specific_rbps) + list(upstream_rbps),
                     gene_sets=['MSigDB_Hallmark_2020'],
                     organism='Human',
                     description='enrichr',
                     outdir='test/enrichr_kegg',
                     cutoff=0.05
                     )
print(enr.results.head(200))

