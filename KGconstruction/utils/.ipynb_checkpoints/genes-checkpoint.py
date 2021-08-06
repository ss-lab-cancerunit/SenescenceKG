import pandas as pd
import numpy as np

class GeneDataManager:
    
    # senescence genes filepath should be a dataframe of genes to include in the graph
    def __init__(self, 
                 genes_filepath: str,
                 drop_dups: bool = True,
                 sep: str = ','):
        
        self.genedata = pd.read_csv(genes_filepath,
                                    sep = sep)
        
        self.genedata.fillna('NA', inplace = True)
        
        if drop_dups:
            self.genedata.drop_duplicates(inplace = True)
        
    @staticmethod
    def mapID(key: str, table: pd.DataFrame, column: str) -> str:
        if key in table.index:
            mapped_id = table.loc[key, column]
        else:
            mapped_id = key
        return mapped_id
    
    # function for converting gene IDs in a list of relationships to IDs of a specific type
    def mapRelationshipIDs(self, relationships: list, keytype: str, column: str) -> list[tuple[str, dict, str]]:
        
        table = self.genedata.drop_duplicates(keytype, keep = 'first').copy()
        table.set_index(keytype, inplace = True)
        
        new_relationships = []
        for head, link, tail in relationships:
            new_head = self.mapID(head, table, column)
            new_tail = self.mapID(tail, table, column)
            new_rel = (new_head, link, new_tail)
            new_relationships.append(new_rel)
            
        return new_relationships
    
    # merge senescence gene data with uniprot ids
    # uniprot id mapping file is found here: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
    # args and kwargs for pd.read_csv() function for the file to merge
    def mergeGenesets(self, 
                      filepath_to_merge: str,
                      geneset_on: str = 'ensembl', 
                      other_on: str = 'ensembl',
                      other_include: list = ['UniprotKB-AC'],
                      other_on_column_delim: str = ';',
                      *args, **kwargs):
        
        merge_data = pd.read_csv(filepath_to_merge, *args, **kwargs)
        
        # clean column to merge on - remove NAs and expand if necessary
        merge_data = merge_data[~merge_data[other_on].isna()]
        if other_on_column_delim:
            merge_data[other_on] = [ens.split(other_on_column_delim) for ens in merge_data[other_on].values]
            merge_data = merge_data.explode(other_on)
        merge_data.drop_duplicates(other_on, keep = 'first', inplace = True)
        
        # subset relevant columns for merging
        cols = other_include + [other_on] if other_on not in other_include else other_include
        to_merge = merge_data[cols]
        
        # merge with senescence genes
        self.genedata = self.genedata.merge(to_merge,
                                            left_on = geneset_on, 
                                            right_on = other_on, 
                                            how = 'left')
    
    # get a dictionary of gene properties from the table
    def getGeneProperties(self, 
                          id_col: str = 'ensembl',
                          ensembl_col: str = 'ensembl',
                          entrez_col: str = 'entrez',
                          symbol_col: str = 'symbol',
                          oncogene_LFC_col: str = 'LFC_Oncogene',
                          replicative_LFC_col: str = 'LFC_Replicative',
                          dnadamage_LFC_col: str = 'LFC_DNAdamage',
                          oncogene_p_col: str = 'adj_p_Oncogene',
                          replicative_p_col: str = 'adj_p_Replicative',
                          dnadamage_p_col: str = 'adj_p_DNAdamage',
                          diffexp_col: str = 'sig_types',
                          source_col: str = 'source',
                          senescence_effect_col: str = 'senescence_effect') -> dict:
        
        properties = {}
        for i, row in self.genedata.iterrows():
            properties[row[id_col]] = {'type': 'Gene',
                                       'properties': {
                                            'ensembl': row[ensembl_col],
                                            'entrez': int(row[entrez_col]) if row[entrez_col] != 'NA' else row[entrez_col],
                                            'symbol': row[symbol_col],
                                            'oncogene_LFC': row[oncogene_LFC_col],
                                            'replicative_LFC': row[replicative_LFC_col],
                                            'DNAdamage_LFC': row[dnadamage_LFC_col],
                                            'oncogene_adj_p': row[oncogene_p_col],
                                            'replicative_adj_p': row[replicative_p_col],
                                            'DNAdamage_adj_p': row[dnadamage_p_col],
                                            'differential_expression': row[diffexp_col], 
                                            'source': row[source_col],
                                            'senescence_effect': row[senescence_effect_col]
                                           }
                                      }
            
        return properties