# An integrative knowledge graph of cellular senescence

### Welcome to my master's thesis! 

This project was undertaken as a member of the [Samarajiwa lab](https://www.samarajiwa-lab.org/) at the Cambridge MRC cancer unit.
Thank you to the entire group, who have been extremely welcoming, supportive, and helpful
during my tenure here. Special thanks to Dr. Samarajiwa himself, who has been an invaluable mentor
throughout the process.

This project explores the biological mechanisms of cellular ageing - also called "senescence" -
using a purpose-built database constructed 
from publicly-available data sets and resources related to senescence. 
We collected gene expression data from 20 different senescence-related experiments 
and used a statistical model to identify genes with senescence-associated expression patterns 
that are consistently reproduced by different experiments.
We used these senescence-associated genes to construct a network database of relationships 
between genes (including physical interactions, gene annotations, etc.). 
Then, using data from the database, we implement machine learning models to 
predict novel biological facts related to cellular senescence. 

### Reproducible research 

For reproducibility purposes, all code and most data files are included here,
with one exception being data files from [Pathway Commons](https://www.pathwaycommons.org/), 
[Gene Ontology](http://geneontology.org/), [ReMap](https://remap.univ-amu.fr/),
and [DrugBank](https://go.drugbank.com/), which are too large to feasibly include in this repository. 

However, these files can be easily obtained from the following direct download links (download will start immediately after opening): 
- [Pathway Commons](https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.txt.gz) 
- [Gene Ontology OWL](http://purl.obolibrary.org/obo/go.owl)
- [Gene Ontology annotations](http://geneontology.org/gene-associations/goa_human.gaf.gz)
- [ReMap](https://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_all_macs2_hg38_v1_0.bed.gz) (large file, ~3.6Gb)
- [DrugBank](https://go.drugbank.com/releases/latest) (downloading requires DrugBank approval)

For similar reasons, we also omit raw sequencing read files (FASTA files etc) for RNA-sequencing data sets, 
though gene-level read counts and DESeq2 results for each treatment/control comparison are included in [`RNAseq/data`](RNAseq/data).

The entire senescence knowledge graph has also been included in json format at
[`KGembedding/data/KnowledgeGraph/SenescenceKG.json`](KGembedding/data/KnowledgeGraph/SenescenceKG.json),
as well as the trained entity embeddings for five knowledge graph embedding models, which are in [`KGembedding/data/models`](KGembedding/data/models), 
stored as pickled python dictionaries. Conda virtual environment .yml files have been added to the three python project directories
[KGanalysis](KGanalysis/envs), [KGembedding](KGembedding/envs), and [KGconstruction](KGconstruction/envs).

## RNA-sequencing meta-analysis

The knowledge graph was constructed with genes identified from a meta-analysis of 33 
senescent/control comparisons taken from 20 experiments. All relevant files
can be found in the [`RNAseq/`](RNAseq) directory. Gene read counts and DESeq2's
estimated log2 fold changes and p-values for individual experiments are included 
in [`RNAseq/data`](RNAseq/data), as well as metadata. 
To identify genes with consistent patterns of differential expression, we use
a linear mixed model to model a gene's estimated log2 fold changes, 
treating study heterogeneity as a random effect:

<img src="https://render.githubusercontent.com/render/math?math=\vec{y}_i = {\bf X}\vec{\beta}_i %2B {\bf Z}\vec{\nu}_i %2B \vec{\epsilon}_i">

\
Where the random effect design matrix **Z** contains dummy variables for experiments
and the fixed effect design matrix **X** contains dummy variables for three subtypes of senescence
(oncogene-induced, replication-induced, and DNA damage-induced). 
Genes with at least one absolute estimated log2 fold change greater than one
and p < 0.05 for the corresponding coefficient (shown in blue below) were added 
to the graph. 

![fig1](figures/fig1b.png)

We also add genes from four public biological databases: [CSgene](http://csgene.bioinfo-minzhao.org/),
[Ageing Atlas](https://ngdc.cncb.ac.cn/aging/index), [CellAge](https://genomics.senescence.info/cells/),
and [InnateDB](https://www.innatedb.com/).

## Constructing the graph

We integrated data from [Pathway Commons](https://www.pathwaycommons.org/), 
[Gene Ontology](http://geneontology.org/), [ReMap](https://remap.univ-amu.fr/),
and [DrugBank](https://go.drugbank.com/) into the knowledge graph, adding relations between
senescence-associated genes (identified from the RNA-sequencing analysis and
the four public gene databases). All data sets are publicly available and can be downloaded from the links provided at the top. 
Code used to parse the data sets and construct the graph can be found in [`KGconstruction/`](KGconstruction). 
An overview of the entities and relations present in the graph is shown below.

![fig2](figures/fig2.png)

## Knowledge graph embedding

Knowledge graph embedding was performed with a custom model ("Hierarchical TransR") 
and four existing models from the [PyKeen](https://github.com/pykeen/pykeen) package: 
TransE, TransR, ConvE, and ComplEx. Hyperparameter searches for all models were 
implemented using [Optuna](https://github.com/optuna/optuna). 
All data sets used and produced in the embedding analysis  
are included in the [`KGembedding/data`](KGembedding/data) directory, 
including the trained embeddings for the five models and the entire knowledge graph's json file.
Training was conducted on NVIDIA TITAN V (12 Gb) and NVIDIA TITAN RTX (24 Gb) GPUs. 
Models were evaluated on the test set using the [adjusted mean rank index](https://arxiv.org/pdf/2002.06914.pdf).
Results are shown below.

![fig3](figures/fig3.png)

## Gene classification

Using the trained embeddings, gene classification models were implemented for six target variables.
Models were fit to three target variables indicating categories of expression levels in the three 
senescence subtypes and three other target variables indicating gene pathway memberships. 
Code for classification model building can be found in the [`KGanalysis/`](KGanalysis) directory.
ROC curves for all models are shown below.

![fig5](figures/fig5.png)


