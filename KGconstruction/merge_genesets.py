import pandas as pd
import numpy as np
import mygene
import re
import argparse
mg = mygene.MyGeneInfo()

# filepaths for database files taken from senescence/immunity gene databases
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cellAge', type = str, help = 'Filepath for semicolon-delimited file downloaded from cellAge')
parser.add_argument('-s', '--csGene', type = str, help = 'Filepath for tab-delimited file downloaded from CS gene')
parser.add_argument('-a', '--ageingAtlas', type = str, help = 'Filepath for comma-delimited file downloaded from ageing atlas')
parser.add_argument('-i', '--innateDB_curated', type = str, help = 'Filepath for tab-delimited file downloaded from innateDB')
parser.add_argument('-g', '--innateDB_GO', type = str, help = 'Filepath for comma-delimited file downloaded from innateDB for genes annotated with GO term: "Innate immunity"')
parser.add_argument('-o', '--outpath', type = str, default = 'senescence_db_genes.csv', help = 'Filepath where the merged gene list should be saved')
args = parser.parse_args()

if __name__ == '__main__':

    # import genes from CellAge database
    cellage_genes = pd.read_csv(args.cellAge, sep = ';')
    cellage_genes['source'] = 'CellAge'
    cellage_genes.rename({'gene_name': 'symbol', 'entrezid': 'entrez'},
                           axis=1, inplace=True)
    
    cellage_genes.drop_duplicates('entrez', keep = 'first', inplace = True)
    
    # import genes from Ageing Atlas database
    ageingatlas_genes = pd.read_csv(args.ageingAtlas)
    ageingatlas_genes.rename({'Symbol': 'symbol', 'Gene_ID': 'entrez',
                              'Species': 'organism'}, 
                             axis = 1, inplace = True)
    ageingatlas_genes = ageingatlas_genes[ageingatlas_genes['organism'] == 'Homo sapiens']
    ageingatlas_genes['organism'] = 'Human'
    ageingatlas_genes['source'] = 'Ageing atlas'
    ageingatlas_genes['senescence_effect'] = 'Unknown'
    
    # import genes from CSgene database
    csgene_genes = pd.read_csv(args.csGene, sep = '\t')
    csgene_genes.rename({'GeneSymb': 'symbol', 'GeneID': 'entrez'}, 
                             axis = 1, inplace = True)
    csgene_genes['organism'] = 'Human'
    csgene_genes['source'] = 'CSgene'
    csgene_genes['senescence_effect'] = 'Unknown'

    csgene_genes.drop_duplicates('entrez', keep = 'first', inplace = True)

    # import genes from innatedb database
    innatedb_genes = pd.read_csv(args.innateDB_curated, sep = '\t')
    
    # filter for human genes
    innatedb_genes = innatedb_genes[innatedb_genes['Species'] == 9606]
    
    # convert concatenated gene names into separate rows
    innatedb_genes['Gene Symbol'] = innatedb_genes['Gene Symbol'].str.split(';')
    innatedb_genes = innatedb_genes.explode('Gene Symbol')
    # convert symbols to uppercase to standardise
    innatedb_genes['Gene Symbol'] = innatedb_genes['Gene Symbol'].str.upper()
    
    # import genes with GO annotation "innate immune response"
    innatedb_genes_GO = pd.read_csv(args.innateDB_GO)
    innatedb_genes_GO = innatedb_genes_GO[innatedb_genes_GO['taxonId'] == 9606]
    innatedb_genes_GO['entrez'] = innatedb_genes_GO['entrez'].str.extract(',*([^,]+)$')

    print('Datasets imported, mapping gene IDs')
    
    # use GO df to get a mapping between symbols and entrez IDs
    GO_gene_entrez_ids = {row['name'].upper(): row['entrez'] for i, row in innatedb_genes_GO.iterrows()}
    
    # get unique genes
    entrez_id_mapping = {}
    unique_genes = pd.unique(innatedb_genes['Gene Symbol'])
    # query genes in chunks of 10
    for i in range(0, len(unique_genes), 10):
        # run query for 10 genes
        max_idx = min(i+10, len(unique_genes))
        genes_to_query = unique_genes[i:max_idx]
        query = 'symbol:' + ' OR symbol:'.join(genes_to_query)
        hits = mg.query(query, species = 9606)['hits']
        # subset hits that have an entrez ID
        all_gene_hits = {hit['symbol'].upper(): hit for hit in hits if 'entrezgene' in hit}
        for gene in genes_to_query:
            # if gene was found in the query, add the entrez ID to dictionary
            if gene in all_gene_hits:
                gene_hit = all_gene_hits[gene]
                entrez_id_mapping[gene] = gene_hit['entrezgene']
            # otherwise, check to see if the gene's ID was provided by the GO data set
            elif gene in GO_gene_entrez_ids:
                entrez_id_mapping[gene] = GO_gene_entrez_ids[gene]
                
    # fill in entrez IDs using the dictionary
    innatedb_genes['entrez'] = innatedb_genes['Gene Symbol'].apply(
        lambda gene: entrez_id_mapping[gene.upper()] if gene.upper() in entrez_id_mapping else np.nan)
    
    # fill in a few missing ones manually
    innatedb_genes.loc[innatedb_genes['Gene Symbol'] == 'DEFB4', 'entrez'] = '1673'
    innatedb_genes.loc[innatedb_genes['Gene Symbol'] == 'IL1F7', 'entrez'] = '27178'
    innatedb_genes.loc[innatedb_genes['Gene Symbol'] == 'IL8RB', 'entrez'] = '3579'
    innatedb_genes.loc[innatedb_genes['Gene Symbol'] == 'ICOSL', 'entrez'] = '23308'
    innatedb_genes.loc[innatedb_genes['Gene Symbol'] == 'NLRP1A', 'entrez'] = '22861'
    innatedb_genes.loc[innatedb_genes['Gene Symbol'] == 'KDM1', 'entrez'] = '23028'
    
    # drop genes that are still NA for entrez ID, convert to integer
    innatedb_genes = innatedb_genes[~innatedb_genes['entrez'].isna()]
    innatedb_genes['entrez'] = innatedb_genes['entrez'].astype(np.int64)
    
    # add/change column names in immune genes data frame to match senescence df
    innatedb_genes['source'] = 'InnateDB'
    innatedb_genes['senescence_effect'] = 'Unknown'
    innatedb_genes['Species'] = 'Human'
    innatedb_genes.rename({'Gene Symbol': 'symbol', 'Species': 'organism'},
                        axis = 1, inplace = True)
    
    innatedb_genes.drop_duplicates('entrez', keep = 'first', inplace = True)
    
    datasets = [cellage_genes, ageingatlas_genes, csgene_genes, innatedb_genes]

    cols = ['symbol', 'entrez', 'organism', 'source', 'senescence_effect']
    for i, dataset in enumerate(datasets):
        if i == 0:
            all_genes = dataset[cols].copy()
        else:
            source = dataset['source'].values[0]
            genes_in_set = all_genes['entrez'].isin(dataset['entrez'])
            all_genes.loc[genes_in_set, 'source'] = all_genes.loc[genes_in_set, 'source'] + f'/{source}'
            new_genes = dataset.loc[~dataset['entrez'].isin(all_genes['entrez']), cols]
            all_genes = pd.concat([all_genes, new_genes], axis = 0)

    # save all genes file to csv
    print('Saving gene info file to: {}'.format(args.outpath))
    all_genes.to_csv(args.outpath, index = False)
        