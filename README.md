# An integrative knowledge graph of cellular senescence

### Welcome to my master's thesis! 

I undertook this project as a member of the [Samarajiwa lab](https://www.samarajiwa-lab.org/) at the Cambridge MRC cancer unit.
Thank you to the entire group, who have been extremely welcoming, supportive, and helpful
during my time here. Special thanks to Dr. Samarajiwa himself, who has been an invaluable mentor
throughout this entire process.

## Non-technical summary

This project explores the biological mechanisms of cellular ageing - also known as "senescence" -
using a purpose-built database, constructed 
from publicly-available data sets collected in experiments studying senescent cells. 
We collected gene expression data from 20 different experiments 
and used a statistical model to identify genes with senescence associations
that were consistently reproduced across different experiments.
We constructed a network database of relationships 
between these genes (including physical interactions, gene annotations, etc.),
hoping to reveal interesting properties of the biological networks that govern the ageing process. 
Using facts from the database, we implemented machine learning models to
predict novel biological facts related to cellular senescence.

## Reproducible research 

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
[`KGanalysis/`](KGanalysis/envs), [`KGembedding/`](KGembedding/envs), and [`KGconstruction/`](KGconstruction/envs).

## Project overview

```{python}
.
├── KGanalysis      # Directory with code for training supervised learning models on gene embeddings
    ├── data        # Directory storing data needed to train classifiers
        ├── embeddings      # Directory storing KG embeddings as pickle files
        ├── classifier      # Directory storing hyperparameter search results
        └── genesets        # Directory storing gene sets used in classification tasks
    ├── envs        # Directory storing conda environment .yml file
    ├── classifier_hp_search.py     # Script for running classifier hyperparameter search 
    └── GeneClassifier.ipynb        # Jupyter notebook for training and visualisation final classification models
├── KGconstruction      # Directory with code for constructing the knowledge graph
    ├── data        # Directory storing data needed to construct the KG
        ├── TFs     # Directory storing gene promoter data
        └── genesets        # Directory storing gene sets used to construct the graph   
    ├── envs        # Directory storing conda environment .yml file
    ├── utils       # Directory with scripts defining classes used to parse database files
    ├── build_graph.py      # Script for building the graph, takes database file names as command line arguments
    ├── merge_genesets.py       # Script for merging gene sets from the four public senescence gene databases
    └── peaks_to_genes.py       # Script for mapping Roadmap promoter regions to genes
├── KGembedding     # Directory with code for knowledge graph embedding models
    ├── data        # Directory storing data needed to train embeddings
        ├── KnowledgeGraph      # Directory storing the knowledge graph as a json file
        ├── hyperparams     # Directory storing embedding model hyperparamter search results
        └── models      # Directory storing final embedding models and train/test sets
    ├── envs        # Directory storing conda environment .yml file
    ├── utils       # Directory storing embedding model code
        ├── embedding.py        # Code for custom embedding model "Hierarchical TransR"
        ├── preprocessing.py        # Code for parsing and cleaning the knowledge graph data
        └── transR.py       # Code for TransR, used in Hierarchical TransR
    ├── EmbeddingVisualisation.ipynb        # Notebook visualising entity embeddings
    ├── LinkPrediction.ipynb        # Notebook predicting new links with embedding models
    ├── ModelEvaluation.ipynb       # Notebook evaluating embedding models using AMRI
    ├── PlotGraphData.ipynb     # Notebook for plotting data on the knowledge graph's composition
    ├── custom_hp_search.py     # Script for hyperparameter search for hierarchical TransR
    ├── fit_models.py       # Script for fitting final embedding models (after hyperparameter searches)
    └── pykeen_hp_search.py     # Script for hyperparameter searches for PyKeen embedding models
├── RNAseq      # Directory with code for RNA-sequencing meta-analysis
    ├── data        # Directory storing experiment read count data and differential expression results
    ├── figures     # Directory storing figures generated by the RNA-sequencing analysis
    ├── metadata        # Directory storing sample and comparison metadata 
    ├── merge_DE.R      # Script for merging all DESeq2 results into a single data set
    ├── merge_counts.R      # Script for merging all sample read counts into a single data set
    ├── meta_analysis.R     # Script for meta-analysis of all 33 comparisons
    ├── plotting.Rmd        # R markdown notebook plotting results
    └── utils.R     # Script with functions used in the analysis
├── figures     # Directory with all figures
├── README.md       # README file
└── report.pdf      # Final report of the entire project
```

## RNA-sequencing meta-analysis

The knowledge graph was constructed with genes identified from a meta-analysis of 33 
senescent/control comparisons taken from 20 experiments. All relevant files
can be found in the [`RNAseq/`](RNAseq) directory. Gene read counts and DESeq2's
estimated log2 fold changes and p-values for individual experiments are included 
in [`RNAseq/data`](RNAseq/data), as well as comparison and sample metadata in [`RNAseq/metadata`](RNAseq/metadata). 
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


