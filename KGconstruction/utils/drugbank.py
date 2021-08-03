import re
import xml.etree.ElementTree as ET
import pandas as pd

class DrugbankParser:
    
    # arguments are filepath to drugbank XML file as well as column names in the parsed data frame
    def __init__(self, 
                 filepath: str,
                 action_col: str = 'actions',
                 alias_col: str = 'drug_aliases',
                 drugname_col: str = 'drug_name',
                 category_col: str = 'drug_category',
                 approval_col: str = 'approval_status',
                 organism_col: str = 'organism',
                 description_col: str = 'description',
                 db_id_col: str = 'drugbank_id',
                 uniprot_id_col: str = 'uniprot_id'):
        
        # open XML file
        self.tree = ET.parse(filepath)
        self.root = self.tree.getroot()
        
        # set instance variables for column names in drug target data frame
        self.action_col = action_col
        self.alias_col = alias_col
        self.drugname_col = drugname_col
        self.category_col = category_col
        self.organism_col = organism_col
        self.approval_col = approval_col
        self.description_col = description_col
        self.db_id_col = db_id_col
        self.uniprot_id_col = uniprot_id_col
        
        # parse XML file to get drug target information
        self.getTargets()

    # function for getting a dataframe of drugs and their targets
    def getTargets(self):
              
        ns = '{http://www.drugbank.ca}'
        protein_rows = []
        for drug in self.root:
            
            # get name and ID
            drugbank_id = drug.findtext(ns + "drugbank-id[@primary='true']")
            name = drug.findtext(ns + 'name').lower()
            
            # get aliases as well
            international_brands = {elem.text.lower() for elem in 
                                    drug.findall('{ns}international-brands/{ns}international-brand/{ns}name'.format(ns = ns))}
            synonyms = {elem.text.lower() for elem in 
                        drug.findall('{ns}synonyms/{ns}synonym'.format(ns=ns))}
            products = {elem.text.lower() for elem in 
                        drug.findall('{ns}products/{ns}product/{ns}name'.format(ns = ns))}
            aliases = international_brands.union(synonyms, products)
            
            drug_categories = [cat.findtext(ns + 'category') for cat in 
                               drug.findall("{ns}categories/{ns}category".format(ns = ns))]
            
            drug_approval_groups = [group.text for group in drug.findall("{ns}groups/{ns}group".format(ns = ns))]
            
            drug_description = drug.findtext(ns + "description")
            
            row = {self.db_id_col: drugbank_id, 
                   self.drugname_col: name, 
                   self.alias_col: '; '.join(aliases), 
                   self.category_col: '; '.join(drug_categories),
                   self.approval_col: '; '.join(drug_approval_groups),
                   self.description_col: drug_description}
            
            # get information on proteins linked to drug
            for category in ['target', 'enzyme', 'carrier', 'transporter']:
                proteins = drug.findall('{ns}{cat}s/{ns}{cat}'.format(ns=ns, cat=category))
                for protein in proteins:
                    row['category'] = category
                    row[self.organism_col] = protein.findtext('{}organism'.format(ns))
                    actions = protein.findall('{ns}actions/{ns}action'.format(ns=ns))
                    row[self.action_col] = ';'.join(action.text for action in actions)
                    uniprot_ids = [polypep.text for polypep in protein.findall(
                        "{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier".format(ns=ns))]            
                    if len(uniprot_ids) != 1:
                        continue
                    row[self.uniprot_id_col] = uniprot_ids[0]
                    protein_rows.append(row)
        
        # put into data frame
        protein_df = pd.DataFrame.from_dict(protein_rows)
        # fill rows with unknown actions
        protein_df.loc[protein_df[self.action_col] == '', self.action_col] = 'unknown'
        # split drug actions into separate rows
        protein_df[self.action_col] = [action.split(';') for action in protein_df[self.action_col].values]
        protein_df = protein_df.explode(self.action_col)
        # drop duplicate rows
        protein_df.drop_duplicates([self.db_id_col, self.uniprot_id_col, self.action_col], inplace = True)
        
        self.targets = protein_df
    
    # get drug-protein relationships
    def getRelationships(self, 
                         uniprot_ids: pd.Series = None) -> list[tuple[str, dict, str]]:

        if uniprot_ids is not None:
            targs = self.targets[self.targets[self.uniprot_id_col].isin(uniprot_ids)].copy()
        else:
            targs = self.targets.copy()
            
        relationships = []
        for i, row in targs.iterrows():
            drug = row[self.db_id_col] 
            gene = row[self.uniprot_id_col]
            target_mech = row[self.action_col].upper().replace(' ', '_') if 'unknown' not in row[self.action_col] else 'TARGETS'
            link = {'type': target_mech,
                    'properties': {
                        'source': 'Drugbank', 
                        }
                    }
            relationship = (drug, link, gene)
            relationships.append(relationship)
            
        return relationships
    
    def getDrugProperties(self, id_col: str = 'drugbank_id') -> dict:
        
        properties = {}
        unique_drugs = self.targets.drop_duplicates([self.db_id_col])
        for i, row in unique_drugs.iterrows():
            properties[row[id_col]] = {'type': 'Drug',
                                       'properties': {
                                            'drugbank_id': row[self.db_id_col],
                                            'organism': row[self.organism_col],
                                            'name': row[self.drugname_col],
                                            'aliases': row[self.alias_col],
                                            'categories': row[self.category_col],
                                            'description': row[self.description_col],
                                            'approval_status': row[self.approval_col]
                                           }
                                      }
            
        return properties         