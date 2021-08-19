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


# In[2]:


sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# # SA609 Patient Analysis

# In[3]:


adata = sc.read("sa609.h5ad", cache=True)
adata.X = csr_matrix(adata.X)


# In[4]:


genes = adata.var_names
gene_list = genes.tolist()


# In[5]:


with open('RBP2GO_RBP_list.csv', newline='', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)
    RBP_list = list(itertools.chain.from_iterable(data))
valid_rbps = set(RBP_list).intersection(set(gene_list))


# In[6]:


sc.pl.umap(adata, color=['clone'])
sc.tl.rank_genes_groups(adata, 'clone', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[7]:


sc.tl.rank_genes_groups(adata, 'therapy', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[8]:


I = pickle.load(open("intron_network.pkl","rb"))


# In[9]:


valid_rbps = set()
for x in I.edges():
    valid_rbps.add(x[0])


# In[10]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subI = I.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subI, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subI, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[11]:


U = pickle.load(open("utr_network.pkl","rb"))


# In[12]:


valid_rbps = set()
for x in U.edges():
    valid_rbps.add(x[0])


# In[13]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subU = U.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subU, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subU, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[14]:


N = pickle.load(open("non_network.pkl","rb"))


# In[15]:


valid_rbps = set()
for x in N.edges():
    valid_rbps.add(x[0])


# In[16]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subN = N.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subN, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subN, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[17]:


S = pickle.load(open("similar_network.pkl","rb"))


# In[18]:


valid_rbps = set()
for x in S.edges():
    valid_rbps.add(x[0])


# In[19]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subS = S.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subS, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subS, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[20]:


rbp_exp = adata.X.todense()[:,adata.var.index.tolist().index("EIF4B")].T.tolist()[0]
target_exp = adata.X.todense()[:,adata.var.index.tolist().index("RBMXL1")].T.tolist()[0]
plt.scatter(target_exp, rbp_exp)


# # SA535 Patient Analysis

# In[21]:


adata = sc.read("sa535.h5ad")
adata.X = csr_matrix(adata.X)


# In[22]:


genes = adata.var_names
gene_list = genes.tolist()


# In[23]:


with open('RBP2GO_RBP_list.csv', newline='', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)
    RBP_list = list(itertools.chain.from_iterable(data))
valid_rbps = set(RBP_list).intersection(set(gene_list))


# In[24]:


sc.pl.umap(adata, color=['clone'])
sc.tl.rank_genes_groups(adata, 'clone', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[25]:


sc.tl.rank_genes_groups(adata, 'therapy', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[26]:


I = pickle.load(open("intron_network.pkl","rb"))


# In[27]:


valid_rbps = set()
for x in I.edges():
    valid_rbps.add(x[0])


# In[28]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subI = I.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subI, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subI, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[29]:


U = pickle.load(open("utr_network.pkl","rb"))


# In[30]:


valid_rbps = set()
for x in U.edges():
    valid_rbps.add(x[0])


# In[31]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subU = U.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subU, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subU, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[32]:


N = pickle.load(open("non_network.pkl","rb"))


# In[33]:


valid_rbps = set()
for x in N.edges():
    valid_rbps.add(x[0])


# In[34]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subN = N.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subN, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subN, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[35]:


S = pickle.load(open("similar_network.pkl","rb"))


# In[36]:


valid_rbps = set()
for x in S.edges():
    valid_rbps.add(x[0])


# In[37]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subS = S.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subS, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subS, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[38]:


rbp_exp = adata.X.todense()[:,adata.var.index.tolist().index("EIF4B")].T.tolist()[0]
target_exp = adata.X.todense()[:,adata.var.index.tolist().index("RBMXL1")].T.tolist()[0]
plt.scatter(target_exp, rbp_exp)


# # SA1035 Patient Analysis

# In[39]:


adata = sc.read("sa1035.h5ad")
adata.X = csr_matrix(adata.X)


# In[40]:


genes = adata.var_names
gene_list = genes.tolist()


# In[41]:


with open('RBP2GO_RBP_list.csv', newline='', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)
    RBP_list = list(itertools.chain.from_iterable(data))
valid_rbps = set(RBP_list).intersection(set(gene_list))


# In[42]:


sc.pl.umap(adata, color=['clone'])
sc.tl.rank_genes_groups(adata, 'clone', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[43]:


sc.tl.rank_genes_groups(adata, 'therapy', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False)


# In[44]:


I = pickle.load(open("intron_network.pkl","rb"))


# In[45]:


valid_rbps = set()
for x in I.edges():
    valid_rbps.add(x[0])


# In[46]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subI = I.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subI, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subI, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[47]:


U = pickle.load(open("utr_network.pkl","rb"))


# In[48]:


valid_rbps = set()
for x in U.edges():
    valid_rbps.add(x[0])


# In[49]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subU = U.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subU, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subU, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[50]:


N = pickle.load(open("non_network.pkl","rb"))


# In[51]:


valid_rbps = set()
for x in N.edges():
    valid_rbps.add(x[0])


# In[52]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subN = N.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subN, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subN, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[53]:


S = pickle.load(open("similar_network.pkl","rb"))


# In[54]:


valid_rbps = set()
for x in S.edges():
    valid_rbps.add(x[0])


# In[55]:


degs = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
for cluster in degs:
    print(cluster)
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
                         gene_sets=['KEGG_2021_Human','MSigDB_Hallmark_2020'],
                         organism='Human',
                         description='rbps',
                         outdir='enrichment',
                         cutoff=0.05
                        )
        for x,y in zip(enr.results["Term"][:50],enr.results["Genes"][:50]):
            if "Ox" in x or "DNA" in x or "cancer" in x:
                print(x,y)
    except Exception as e:
        continue
    subS = S.edge_subgraph(valid_edges)
    pos = nx.nx_agraph.graphviz_layout(subS, prog="neato",args="-Goverlap=scale -Elen=5 -Eweight=0.2")
    fig = plt.figure(figsize=(16,8))
    ax = plt.subplot(1,2,1)
    nx.draw(subS, pos, with_labels=True, ax=ax)
    ax = plt.subplot(1,2,2)
    sc.pl.umap(adata, ax=ax,show=False)
    cluster_sub_umap = adata[adata.obs["therapy"] == cluster]
    sc.pl.umap(cluster_sub_umap,color="therapy",ax=ax,show=False)
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# In[56]:


rbp_exp = adata.X.todense()[:,adata.var.index.tolist().index("EIF4B")].T.tolist()[0]
target_exp = adata.X.todense()[:,adata.var.index.tolist().index("RBMXL1")].T.tolist()[0]
plt.scatter(target_exp, rbp_exp)

