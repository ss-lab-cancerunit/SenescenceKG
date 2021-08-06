import pandas as pd
import re

class PathwayCommonsParser:
    
    # initialise by specifying filepaths and column names of the interaction data set to be read-in
    # interactions file is taken from pathway commons, while names/hierarchy files are taken from reactome
    def __init__(self, 
                 interactions_filepath: str, 
                 filesep: str = '\t',
                 protein1_col: str = 'PARTICIPANT_A', 
                 protein2_col: str = 'PARTICIPANT_B', 
                 interaction_col: str = 'INTERACTION_TYPE',
                 interactions_to_filter: list = ['Reference'], 
                 sources_to_keep: list = ['Reactome', 'BioGRID', 'KEGG'],
                 pathway_col: str = 'PATHWAY_NAMES',
                 source_col: str = 'INTERACTION_DATA_SOURCE', 
                 species: str = 'Homo sapiens'):
        
        # import the three data sets provided as arguments
        self.interactions = pd.read_csv(interactions_filepath, sep = filesep)
        
        # set values of column names for interaction data
        self.protein1_col = protein1_col
        self.protein2_col = protein2_col
        self.interaction_col = interaction_col
        self.sources_to_keep = sources_to_keep
        
        # interactions to filter should be a list or set 
        self.interactions_to_filter = interactions_to_filter
        self.pathway_col = pathway_col
        self.source_col = source_col
        
        self.cleanData()
        
    # clean pathway data - keep only pathways from the specified database
    def cleanData(self, pathway_sep: str = ';'):
        
        # keep only pathways from the specified sources
        interactions = self.interactions[~self.interactions[self.source_col].isna()].copy()
        interactions[self.source_col] = [';'.join([src for src in sources.split(';') if src in self.sources_to_keep])
                                         for sources in interactions[self.source_col].values]
        interactions = interactions[interactions[self.source_col] != '']
        
         # annotate NAs in pathway column as "physical interaction"
        interactions.loc[interactions[self.pathway_col].isna(), self.interaction_col] = 'PHYSICAL_INTERACTION'
        interactions.loc[interactions[self.pathway_col].isna(), self.pathway_col] = 'NA'
        
        # aggregate duplicate entries from different databases
        interactions = interactions.groupby([self.protein1_col, self.interaction_col, 
                                             self.protein2_col, self.pathway_col]).agg(';'.join)
        
        interactions.reset_index(drop = False, inplace = True)
        
        # remove interactions specified at initialisation
        interactions = interactions[~interactions[self.interaction_col].str.contains('|'.join(self.interactions_to_filter))]
        interactions = interactions.iloc[:-1, :]
        
        # make uppercase to standardise
        interactions[self.protein1_col] = interactions[self.protein1_col].str.upper()
        interactions[self.protein2_col] = interactions[self.protein2_col].str.upper()
        
        # convert rows with multiple pathways into multiple rows
        interactions[self.pathway_col] = interactions[self.pathway_col].str.split('(?<=[^\s]);(?=[^\s])')
        interactions = interactions.explode(self.pathway_col)
        
        # remove NA rows
        interactions = interactions[~(interactions[self.protein1_col].isna() | interactions[self.protein2_col].isna())]
        # convert dashes to underscores to standardise
        interactions[self.interaction_col] = [intr.replace('-', '_') for intr in interactions[self.interaction_col]]
        
        self.interactions_cleaned = interactions
        
    # returns all pathway relationships (or only those involving a specified subset of proteins)
    def getGeneRelationships(self, gene_symbols: pd.Series = None) -> list[tuple[str, dict, str]]:
        
        # if genes of interest are specified, subset pathways
        if gene_symbols is not None:
            pathways_of_interest = (self.interactions_cleaned[self.protein1_col].isin(gene_symbols) & 
                                    self.interactions_cleaned[self.protein2_col].isin(gene_symbols))
            relevant_pathways = self.interactions_cleaned[pathways_of_interest].copy()
        else:
            relevant_pathways = self.interactions_cleaned

        # make list of relationships for the graph
        relationships = []
        for i, row in relevant_pathways.iterrows():
            
            head = row[self.protein1_col]
            link = {'type': row[self.interaction_col].upper(),
                    'properties': { 
                        'source': '; '.join([source for source in row[self.source_col].split(';') if source in self.sources_to_keep]) 
                        }
                    }
            
            if row[self.pathway_col] != 'NA':
                link['properties']['pathway_name'] = row[self.pathway_col]

            tail = row[self.protein2_col]
            rel = (head, link, tail)
            relationships.append(rel)
            
            # for physical interactions, add the symmetric relation to the list
            if row[self.interaction_col].upper() == 'PHYSICAL_INTERACTION':
                relationships.append((tail, link, head))
            
        return relationships