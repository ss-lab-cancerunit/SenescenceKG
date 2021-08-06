import numpy as np
import pandas as pd

import argparse

from utils.drugbank import DrugbankParser
from utils.GO import GOParser
from utils.graph import GraphAPI
from utils.pathways import PathwayCommonsParser
from utils.remap import RemapParser
from utils.genes import GeneDataManager

# add command line arguments for data set files
parser = argparse.ArgumentParser()
parser.add_argument('--GO_owl', type = str, help = 'Filepath to Gene Ontology .owl file')
parser.add_argument('--GO_ann', type = str, help = 'Filepath to Gene Ontology annotation .gaf file')
parser.add_argument('--drugbank', type = str, help = 'Filepath to DrugBank .xml file')
parser.add_argument('--pathway_commons', type = str, help = 'Filepath to Pathway Commons .txt file')
parser.add_argument('--remap', type = str, help = 'Filepath to ReMap .bed file')
parser.add_argument('--uniprot_idmapping', type = str, help = 'Filepath to Uniprot .tab file for id mapping')
parser.add_argument('--password', type = str, help = 'Password for the neo4j graph database')
args = parser.parse_args()

# initialise class instances for data parsers
GeneOntology = GOParser(args.GO_owl, args.GO_ann)
DrugBank = DrugbankParser(args.drugbank)
Pathways = PathwayCommonsParser(args.pathway_commons)
Remap = RemapParser(args.remap, 'data/TFs/gene_promoters.tsv')

# import senescence gene data
senescence_genes = GeneDataManager('data/genesets/all_senescence_genes.csv')
senescence_genes.mergeGenesets(args.uniprot_idmapping,
                               geneset_on = 'ensembl', 
                               other_on = 'ensembl',
                               sep = '\t', header = None, 
                               names = ['UniprotKB-AC', 'UniprotKB-ID', 'entrez', 
                                        'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 
                                        'UniRef90', 'UniRef50', 'UniParc', 
                                        'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 
                                        'PubMed', 'EMBL', 'EMBL-CDS', 'ensembl',
                                        'ensembl_TRS', 'ensembl_PRO', 'Additional_PubMed'],
                               low_memory = False)

# get gene symbols and uniprot IDs for senescence genes
senescence_ensembl_ids = senescence_genes.genedata['ensembl'].dropna()
senescence_gene_symbols = senescence_genes.genedata['symbol'].dropna()
senescence_uniprot_ids = senescence_genes.genedata['UniprotKB-AC'].dropna()

# get gene relationships with GO terms
GO_gene_rels = GeneOntology.getGeneRelationships(senescence_gene_symbols)
GO_gene_rels = senescence_genes.mapRelationshipIDs(GO_gene_rels, 'symbol', 'ensembl')
GO_terms_included = set([term for gene, rel, term in GO_gene_rels])

# get relationships between GO terms
GO_term_rels = GeneOntology.getTermRelationships()

# filter for terms that are in the graph
GO_term_rels = [(h, l, t) for h, l, t in GO_term_rels 
                if h in GO_terms_included and t in GO_terms_included]

# same for drugbank
DB_gene_rels = DrugBank.getRelationships(senescence_uniprot_ids)
DB_gene_rels = senescence_genes.mapRelationshipIDs(DB_gene_rels, 'UniprotKB-AC', 'ensembl')

# same for pathway commons
PC_gene_rels = Pathways.getGeneRelationships(senescence_gene_symbols)
PC_gene_rels = senescence_genes.mapRelationshipIDs(PC_gene_rels, 'symbol', 'ensembl')

# same for remap
remap_gene_rels = Remap.getRelationships(senescence_gene_symbols, senescence_ensembl_ids)
remap_gene_rels = senescence_genes.mapRelationshipIDs(remap_gene_rels, 'symbol', 'ensembl')

# combine all relationships into a single list
all_relationships = GO_gene_rels + GO_term_rels + DB_gene_rels + PC_gene_rels + remap_gene_rels

# get properties
gene_properties = senescence_genes.getGeneProperties()
GO_properties = GeneOntology.getTermProperties()
drug_properties = DrugBank.getDrugProperties()

# combine all properties into single list
all_properties = {**gene_properties, **GO_properties, **drug_properties}

# get genes not in relationships
genes_in_rels = set([h for h, l, t in all_relationships if all_properties[h]['type'] == 'Gene'] + 
                    [t for h, l, t in all_relationships if all_properties[t]['type'] == 'Gene'])
genes_not_in_rels = [ent for ent in all_properties if all_properties[ent]['type'] == 'Gene' and ent not in genes_in_rels]

if __name__ == '__main__':
    
    print('\nProperties and relationships generated, building graph...\n')
    
    uri = 'bolt://localhost:7687'
    user = 'neo4j'
    password = args.password
    graph_db = GraphAPI(uri, user, password)
    graph_db.writeGraph(all_relationships, 
                        all_properties,
                        genes_not_in_rels)
    graph_db.close()