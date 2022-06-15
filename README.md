# Modeling RNA Binding Protein Regulatory Interactions and Expression Networks within Ovarian Cancer and Breast Cancer Single Cell RNA Sequencing (scRNA) Datasets

**_Abstract_**: High Grade Serous Ovarian Cancer (HGSOC) is one of the deadliest forms of ovarian cancer, leading as the 8th most common cause of cancer-related deaths among women worldwide. Genomically, HGSOC is characterized by p53 mutations and instability within pathways pertaining to DNA repair. Triple Negative Breast Cancer (TNBC) is representative of 10%-15% of all diagnosed breast cancer patients and is known for its manifestation of chemotherapeutic cisplatin resistance. RNA Binding Proteins (RBP) are proteins that modulate several cellular pathways and structures, via post-transcriptional and post-translational regulatory interactions. Dysregulated RBPs have shown to promote the metastasis and growth of both ovarian and breast cancers, thus suggesting the need for greater knowledge of RBP-target interactions and expression models within HGSOC and TNBC. Furthermore, RBP regulation has shown great involvement in cisplatin resistance within TNBC, thus reasoning the exploration of RBP-targets in cisplatin-treated samples. scRNA sequencing is a genomic approach involved in producing quantitative sequence profiles from mRNA transcripts. ScanPy, a Python-based package for analyzing Single-Cell RNA (scRNA) Seq, was utilized for analysis of ©MSKCC Isabel acquired patient data with networks developed from MotifMap-RNA acquired RBP datasets. Gene enrichment analysis was also done using GSEApy, to find probabilistic cellular pathways for all genes of interest in HGSOC and TNBC. Overall, 81 novel HGSOC-related RBP target sets were identified from the 6 RBPs of interest discovered within all 41 HGSOC patients. Of the 6 RPBs, 3 proteins were also derived from the IGF2 family which have shown great contribution in increasing the cellular incidence of ovarian cancer. Additionally, identified pathways from TNBC were cross-checked using previous research, to produce an array of targets for every RBP of interest within all 3 TNBC patients.

**Necessary Dependencies Required to Install**
1. Python (Version 3.9.5)
2. ScanPy (Version 1.8.1)
3. GSEApy (Version 0.10.5)
4. NetworkX (Versions 2.6.2)
5. SciPy (Version 1.7.0)
6. NumPy (Version 1.20.2)
7. Pandas (Version 1.2.5)
8. PyGraphviz (Version 1.7)
