import pandas as pd
import numpy as np
from numba import jit

class RemapParser:
    
    # both the remap and promoter files should be tab-separated
    def __init__(self, remap_filepath: str, promoter_filepath: str,
                 sep_remap: str = '\t', sep_promoters: str = '\t',
                 score_threshold: int = 2):
        
        # import remap data
        remap = pd.read_csv(remap_filepath, sep = sep_remap,
                            names = ['Chr', 'Start', 'End', 'TF:Cell line', 'Score', 
                                     'Strand', 'thickStart', 'thickEnd', 'RGB'])
        
        # filter by score threshold 
        remap = remap[remap['Score'] >= score_threshold]
        
        # split the TF:cell line column
        remap[['TF', 'Cell line']] = pd.DataFrame([TFCL.split(':') for TFCL in remap['TF:Cell line'].values])
        remap.drop('TF:Cell line', axis = 1, inplace = True)
        remap.drop_duplicates(['TF', 'Chr', 'Start', 'End'], inplace = True)

        # import promoter file
        promoters = pd.read_csv(promoter_filepath, sep = sep_promoters, header = 0)
        promoters.rename({'start_position': 'Start', 'end_position': 'End'}, axis = 1, inplace = True)
        promoters.drop_duplicates(['feature', 'Start', 'End', 'feature_strand'], inplace = True)
        
        self.TFBS = remap.sort_values(['Start', 'End'])
        self.TFBS_by_TF = {TF: sites for TF, sites in self.TFBS.groupby('TF')}
        self.promoters = promoters.sort_values(['Start', 'End'])
        self.promoters_by_chr = {Chrom: promoters for Chrom, promoters in promoters.groupby('chr')}
    
    # function for masking elements in a data frame that are within a specified start, end range
    @staticmethod
    @jit(nopython = True)
    def inRegion(start_sites: np.ndarray, end_sites: np.ndarray, start: int, end: int) -> np.ndarray:
        return (((start_sites > start) & (start_sites < end)) | ((end_sites > start) & (end_sites < end)) | 
                ((start_sites > start) & (end_sites < end)) | ((start_sites < start) & (end_sites > end)))
    
    @staticmethod
    def getSiteTargets(start: int, end: int, promoters: pd.DataFrame) -> list:
        
        overlapping = RemapParser.inRegion(promoters['Start'].to_numpy(), 
                                           promoters['End'].to_numpy(),
                                           start, end)
        
        return promoters.loc[overlapping, 'feature'].to_list()
        
    # function for generating a list of TF binding relationships
    def getRelationships(self, TF_symbols: pd.Series = None, 
                         target_ensembl_ids: pd.Series = None) -> list[tuple[str, dict, str]]:
        
        if TF_symbols is not None:
            genes = [symbol for symbol in TF_symbols if symbol in self.TFBS_by_TF]
        else:
            genes = list(self.TFBS_by_TF.keys())
            
        if target_ensembl_ids is not None:
            promoters = self.promoters[self.promoters['feature'].isin(target_ensembl_ids)]
            promoters_by_chr = {Chrom: promoters for Chrom, promoters in promoters.groupby('chr')}
        else:
            promoters_by_chr = self.promoters_by_chr
        
        # functions for aggregating data from multiple binding sites for the same target
        target_agg_funcs = {'Cell line': set, 'Score': max}
        
        # loop through TFs and make a relationship for each target
        relationships = []
        for TF in genes:
            TFBS = self.TFBS_by_TF[TF]
            # get targets of each TFBS for the current TF
            TFBS['targets']  = TFBS.apply(
                lambda site: self.getSiteTargets(site['Start'], site['End'], 
                                                 promoters_by_chr[site['Chr']]),
                axis = 1
            )
            # remove TFBS with no targets
            TFBS = TFBS[TFBS['targets'].map(lambda targs: len(targs) > 0)]
            
            # only proceed if any targets were found
            if len(TFBS) > 0:
                
                # make a new row for each target
                TFBS = TFBS.explode('targets')
                # aggregate information by target
                TFBS = TFBS[['Cell line', 'Score', 'targets']].groupby('targets').aggregate(target_agg_funcs)
                TFBS.reset_index(inplace = True)
                
                # make a relationship for each TF/target combination
                TFrelationships = TFBS.apply(
                    lambda row: (TF, 
                                 {'type': 'TF_BINDING',
                                  'properties': {
                                      'source': 'ReMap',
                                      'cell_line': '; '.join(row['Cell line']),
                                      'max_score': row['Score']
                                      }
                                 }, 
                                 row['targets']), 
                    axis = 1
                ).to_list()
                
                relationships.extend(TFrelationships)           
            
        return relationships