import owlready2 as OR
import re
import pandas as pd

class GOParser:
    
    # annotation should be in GAF format
    def __init__(self, 
                 owl_filepath: str,
                 annotation_filepath: str, 
                 sep: str = '\t'):
        
        self.GO_annotation = pd.read_csv(annotation_filepath, sep = sep, 
                                         skiprows = 41, low_memory = False, header = None,
                                         names = ['IDtype', 'UniprotID', 'Symbol', 
                                                  'Relationship', 'GO_ID', 
                                                  'DB:Reference', 'Evidence'],
                                         usecols = range(7))
        
        self.go_ontology = OR.get_ontology('file://' + owl_filepath).load()
        
        self.GO_annotation['GO_ID'] = self.GO_annotation['GO_ID'].str.replace(':', '_')
        
        self.getParents()
        self.getDescriptions()
        
         # set function for displaying a GO object
        OR.set_render_func(GOParser.render)

    # make a function that combines GO name and label into single label
    @staticmethod
    def render(entity: str) -> str:
        label = entity.label.first()
        if label:  
            return '{}: {}'.format(entity.name, label)
        else:
            return entity.name
    
    # function for generating a dictionary for all child-parent relationships in the ontology
    def getParents(self):
                
        # fill in dictionary for parents terms of each term
        parents = {}
        
        # loop through terms and add parents to the dictionary
        for term in self.go_ontology.classes():
            terms =  [t for t in list(term.is_a) if 
                      isinstance(t, OR.class_construct.Restriction) 
                      or t._name != 'Thing']
            if terms:
                parents[term.name] = terms
                
        self.parents = parents
        
    def getDescriptions(self):
                
        descriptions = {term.name: '; '.join(term.label)
                        for term in self.go_ontology.classes()
                        if isinstance(term, OR.entity.ThingClass) 
                        and term._name != 'Thing' and term.label}
        
        self.descriptions = descriptions
        
        
    # get list of relationships between GO terms
    def getTermRelationships(self) -> list[tuple[str, dict, str]]:
        
        relationships = []
        # iterate through parent dictionary, save relationships to a list
        for child, parents in self.parents.items():
            for parent in parents:
                # if the parent is just an ontology, save relation as 'is_a'
                if isinstance(parent, OR.entity.ThingClass):
                    link = {'type': 'IS_A', 'properties': {'source': 'Gene ontology'}}
                    relationship = (child, link, parent.name)
                # if the parent object comes with some specific relationship (restriction), save as such
                elif isinstance(parent, OR.class_construct.Restriction):
                    link_type = parent.property.label[0].replace(' ', '_').upper()
                    link = {'type': link_type, 'properties': {'source': 'Gene ontology'}}
                    parent_term = re.sub(': .+$', '', str(parent.value))
                    relationship = (child, link, parent_term)
                # add relationships to the relationship dict
                relationships.append(relationship)
                
        return relationships 
    
    # get list of relationships between genes and GO terms
    def getGeneRelationships(self, gene_symbols: pd.Series = None) -> list[tuple[str, dict, str]]:
                
        if gene_symbols is not None:
            genes = self.GO_annotation[self.GO_annotation['Symbol'].isin(gene_symbols)].copy()
        else:
            genes = self.GO_annotation.copy()
            
        genes = genes.drop_duplicates(['GO_ID', 'Symbol', 'Relationship'])
            
        relationships = []
        for i, row in genes.iterrows():
            if row['GO_ID'] in self.descriptions:
                gene = row['Symbol']
                link = {'type': row['Relationship'].upper().replace('|', '_'), 
                        'properties': {
                            'source': 'Gene ontology', 
                            }
                        }
                relationship = (gene, link, row['GO_ID'])
                relationships.append(relationship)
                
        return relationships
        
    # function for getting the correct GO classification for a given term (from a set of categories)
    @staticmethod
    def getTermClass(term, parent_dict, categories):
        
        current_term = term
        # move up the tree until you hit one of the specified categories
        while current_term not in categories:
            current_term = re.sub(': .+$', '', current_term)
            current_term = str(parent_dict[current_term][0])
            
        return current_term

    # get the classes for each term 
    def getTermProperties(self, 
                          categories = ['GO_0008150: biological_process', 
                                        'GO_0005575: cellular_component', 
                                        'GO_0003674: molecular_function']):
                
        properties = {go_id: {'type': 'GeneOntologyTerm', 
                              'properties': {
                                 'id': go_id,
                                 'description': self.descriptions[go_id],
                                 'class': self.getTermClass(go_id + ': ' + self.descriptions[go_id], 
                                                            self.parents, categories)
                                 }
                             }
                      for go_id in self.descriptions.keys() 
                      if go_id + ': ' + self.descriptions[go_id] in categories 
                      or go_id in self.parents
                     }
        
        return properties