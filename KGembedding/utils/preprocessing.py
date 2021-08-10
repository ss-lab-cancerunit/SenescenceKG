import numpy as np
import pandas as pd
import json
from collections import OrderedDict
from sklearn.model_selection import train_test_split

class GraphParser:
    
    def __init__(self, 
                 json_filepath: str):
        
        # import json file, convert to list of dictionaries
        with open(json_filepath, 'rb') as f:
            jsonlines = [line.decode('utf-8') for line in f.readlines()]
            
        self.jsondata = [json.loads(entry) for entry in jsonlines]
        
        ## extract nodes, relations, and node properties from the json file
        
        # initialise dictionaries/lists for storing graph data
        self.node_entries = {}
        self.node_id2lab = {}
        self.node_lab2id = {}
        self.node_id2idx = {}
        self.relation_entries = {}
        self.general_relations = []
        self.pathway_specific_relations = []
        self.pathway_general_relations = []
        
        # extract data from the json and put into dictionaries
        for entry in self.jsondata:
            if entry['type'] == 'node':
                entry['label'] = '; '.join(entry['labels'])
                node_id = int(entry['id'])
                self.node_entries[node_id] = entry  
                if entry['labels'] == ['Gene']:
                    self.node_id2lab[node_id] = entry['properties']['ensembl']
                    self.node_lab2id[entry['properties']['ensembl']] = node_id
                elif entry['labels'] == ['GeneOntologyTerm']:
                    self.node_id2lab[node_id] = entry['properties']['id']
                    self.node_lab2id[entry['properties']['id']] = node_id
                elif entry['labels'] == ['Drug']:
                    self.node_id2lab[node_id] = entry['properties']['name']
                    self.node_lab2id[entry['properties']['name']] = node_id
            elif entry['type'] == 'relationship':
                self.relation_entries[int(entry['id'])] = entry
                if 'pathway_name' in entry['properties']:
                    self.pathway_general_relations.append(entry['label'])
                    self.pathway_specific_relations.append(entry['properties']['pathway_name'])
                self.general_relations.append(entry['label'])
        
        self.pathway_general_relations = np.unique(self.pathway_general_relations)
        self.pathway_specific_relations = np.unique(self.pathway_specific_relations)
        self.general_relations = np.unique(self.general_relations)
                    
        self.nodes = np.unique(self.node_entries.keys())
        
        # make ID dictionarys for general relations
        self.gen_rel_id2lab = {}
        self.gen_rel_lab2id = {}
        for i, lab in enumerate(self.general_relations):
            self.gen_rel_id2lab[i] = lab
            self.gen_rel_lab2id[lab] = i
            
        # get specific relations
        self.specific_relations = np.concatenate([self.general_relations[~np.isin(self.general_relations, self.pathway_general_relations)], 
                                                  self.pathway_specific_relations])
            
        # make ID dictionarys for specific relations
        self.spec_rel_id2lab = {}
        self.spec_rel_lab2id = {}
        for i, lab in enumerate(self.specific_relations):
            self.spec_rel_id2lab[i] = lab
            self.spec_rel_lab2id[lab] = i
        
        # loop over relations and add to a data frame
        relation_data = {'head': [],
                         'head_lab': [],
                         'head_type': [],
                         'general_relation': [],
                         'general_relation_lab': [],
                         'specific_relation': [],
                         'specific_relation_lab': [],
                         'tail': [],
                         'tail_lab': [],
                         'tail_type': []}
        
        for rel in self.relation_entries.values():
            head = int(rel['start']['id'])
            tail = int(rel['end']['id'])
            head_type = self.node_entries[head]['label']
            tail_type = self.node_entries[tail]['label']
            relation_data['head'].append(head)
            relation_data['head_lab'].append(self.node_id2lab[head])
            relation_data['head_type'].append(head_type)
            relation_data['tail'].append(tail)
            relation_data['tail_lab'].append(self.node_id2lab[tail])
            relation_data['tail_type'].append(tail_type)
            gen_rel = rel['label']
            relation_data['general_relation'].append(self.gen_rel_lab2id[gen_rel])
            relation_data['general_relation_lab'].append(gen_rel)
            # if pathway is specified, add to specific relation
            if 'pathway_name' in rel['properties']:
                spec_rel = rel['properties']['pathway_name']
            # otherwise, specific relation is the general relation
            else:
                spec_rel = gen_rel
            relation_data['specific_relation'].append(self.spec_rel_lab2id[spec_rel])
            relation_data['specific_relation_lab'].append(spec_rel)

        self.graph_data = pd.DataFrame(relation_data)
        
    # clean graph data 
    def processGraph(self,
                     nodes_to_filter: list = ['Drug'],
                     min_num_general: int = 500,
                     min_num_specific: int = 100,
                     state: int = 123):
        
        # filter node types specified
        if nodes_to_filter:
            ids_to_filter = [node_id for node_id in self.node_entries if 
                             any([lab in nodes_to_filter for lab in self.node_entries[node_id]['labels']])]
            
        
        filtered_mask = self.graph_data['head'].isin(ids_to_filter) | self.graph_data['tail'].isin(ids_to_filter)
        graphdata_filtered = self.graph_data[~filtered_mask]
        
        # filter low-frequency relationships according to specified thresholds
        general_relation_counts =  graphdata_filtered['general_relation'].value_counts()
        self.embedding_gen_rels = set(general_relation_counts[general_relation_counts > min_num_general].index)
        specific_relation_counts =  graphdata_filtered['specific_relation'].value_counts()
        self.embedding_spec_rels = set(specific_relation_counts[specific_relation_counts > min_num_specific].index)
        
        rows_to_keep = (graphdata_filtered['general_relation'].isin(self.embedding_gen_rels) 
                        & graphdata_filtered['specific_relation'].isin(self.embedding_spec_rels))
        self.embedding_data = graphdata_filtered[rows_to_keep].copy()
        
        # drop duplicates in case there are any
        self.embedding_data.drop_duplicates(inplace = True)
        
    # function for resolving leaked triples in a train/test split
    @staticmethod
    def fixLeakedTriples(train: pd.DataFrame,
                         test: pd.DataFrame,
                         symmetric_relations: list,
                         relation_column: str,
                         state: int):
        
        # isolate symmetric relations
        train_symmetric_relations = train[train[relation_column].isin(symmetric_relations)]
        test_symmetric_relations = test[test[relation_column].isin(symmetric_relations)]
        
        # get indices of each symmetric relation in train/test sets
        train_triplets_to_idx = OrderedDict({(t, r, h): i for i, (h, r, t) in 
                                             train_symmetric_relations[['head_lab', relation_column, 'tail_lab']].iterrows()
                                            })
        test_triplets_to_idx = OrderedDict({(h, r, t): i for i, (h, r, t) in 
                                            test_symmetric_relations[['head_lab', relation_column, 'tail_lab']].iterrows()
                                           })
        
        # get leaked triplets
        leaked_triplets = [triplet for triplet in test_triplets_to_idx.keys() if triplet in train_triplets_to_idx.keys()]
        # get indices of leaked triplets
        leaked_triplet_idx = {(h, r, t): (train_triplets_to_idx[(h, r, t)], 
                                          test_triplets_to_idx[(h, r, t)])
                              for (h, r, t) in leaked_triplets}
        
        train_to_test = False
        train_new_facts = []
        train_drop_idx = []
        test_new_facts = []
        test_drop_idx = []
        # iterate through leaked triplets and alternate between moving triplets to test and train sets
        for (h, r, t) in leaked_triplets:
            train_idx, test_idx = leaked_triplet_idx[(h, r, t)]
            if train_to_test == True:
                fact = train.loc[train_idx].copy()
                test_new_facts.append(fact)
                train_drop_idx.append(train_idx)
            else:
                fact = test.loc[test_idx].copy()
                train_new_facts.append(fact)
                test_drop_idx.append(test_idx)
            train_to_test = not train_to_test
            
        train = train.drop(train_drop_idx, axis = 0)
        test = test.drop(test_drop_idx, axis = 0)
        train = pd.concat([train, pd.DataFrame(train_new_facts)])
        test = pd.concat([test, pd.DataFrame(test_new_facts)])
        
        train = train.sample(frac = 1, random_state = state)
        test = test.sample(frac = 1, random_state = state)
                
        return train, test
        
    # train/test/validation split
    def split(self, 
              facts: pd.DataFrame,
              test_size: float,
              validation_size: float,
              stratify_column: str = 'specific_relation_lab',
              symmetric_relations: list = ['PHYSICAL_INTERACTION'],
              symmetric_relation_column: str = 'general_relation_lab',
              state: int = 123):
        
         # split data into train/test
        train_set, test_set, _, _ = train_test_split(facts, 
                                                     facts[stratify_column],
                                                     test_size = test_size,
                                                     random_state = state)
        
        train_set, test_set = GraphParser.fixLeakedTriples(train_set, test_set, 
                                                           symmetric_relations = symmetric_relations,
                                                           relation_column = symmetric_relation_column,
                                                           state = state)
        
        # do train/validation split
        train_set, validation_set, _, _ = train_test_split(train_set, 
                                                           train_set[stratify_column],
                                                           test_size = validation_size / (1 - test_size),
                                                           random_state = state) 
        train_set, validation_set = GraphParser.fixLeakedTriples(train_set, validation_set, 
                                                                 symmetric_relations = symmetric_relations,
                                                                 relation_column = symmetric_relation_column,
                                                                 state = state)
        
        # remove entities from the validation set that do not appear in the training set
        train_entities = pd.unique(train_set['head'].append(train_set['tail']))
        validation_set = validation_set[validation_set['head'].isin(train_entities) & validation_set['tail'].isin(train_entities)]
        
        # combine train and validation sets to form the final train set
        final_train_set = pd.concat([train_set, validation_set], axis = 0)
        
        # remove entities from the test set that are not present in the final train set
        final_train_entities = pd.unique(final_train_set['head'].append(final_train_set['tail']))
        test_set = test_set[test_set['head'].isin(final_train_entities) & test_set['tail'].isin(final_train_entities)]

        # get normalised indices for all entities/relations in all data sets
        train_set, final_train_set, test_set, validation_set = self.getEmbeddingIndices(train_set, final_train_set, test_set, validation_set)
        
        self.train_set = train_set
        self.final_train_set = final_train_set
        self.validation_set = validation_set
        self.test_set = test_set
        self.embedding_facts = pd.concat([self.train_set, self.test_set, self.validation_set], axis = 0)
        self.embedding_entities = np.sort(pd.unique(self.embedding_facts['head'].append(self.embedding_facts['tail'])))
        self.num_entities = len(self.embedding_entities)
        
        # use train/test/validation sets to get positive and negative facts in the graph
        self.getFacts()
        
    # function for generating normalised indices over a set of fact data frames (train/test/validation)
    def getEmbeddingIndices(self, *dfs):
        
        all_facts = pd.concat(dfs, axis = 0)
        
        # map entity and relation labels to indices
        self.embedding_node_id2idx = {node: idx for idx, node in enumerate(np.sort(np.unique(all_facts['head'].append(all_facts['tail']))))}
        self.embedding_node_idx2id = {idx: node for node, idx in self.embedding_node_id2idx.items()}
        self.embedding_node_lab2idx = {self.node_id2lab[node]: self.embedding_node_id2idx[node] for node in self.embedding_node_id2idx.keys()}
        self.embedding_node_idx2lab = {idx: lab for lab, idx in self.embedding_node_lab2idx.items()}
        self.embedding_spec_rel_id2idx = {rel: i for i, rel in enumerate(np.sort(np.unique(all_facts['specific_relation'])))}
        self.embedding_gen_rel_id2idx = {rel: i for i, rel in enumerate(np.sort(np.unique(all_facts['general_relation'])))}
        self.embedding_spec_rel_lab2idx = {self.spec_rel_id2lab[rel]: self.embedding_spec_rel_id2idx[rel] 
                                           for rel in self.embedding_spec_rel_id2idx.keys()}
        self.embedding_gen_rel_lab2idx = {self.gen_rel_id2lab[rel]: self.embedding_gen_rel_id2idx[rel]
                                          for rel in self.embedding_gen_rel_id2idx.keys()}
        
        # change entries in the data frame from labels to indices
        for df in dfs:
            df['head_id'] = df['head'].copy()
            df['head'] = df['head'].replace(self.embedding_node_id2idx)
            df['tail_id'] = df['tail'].copy()
            df['tail'] = df['tail'].replace(self.embedding_node_id2idx)
            df['specific_relation'] = df['specific_relation'].replace(self.embedding_spec_rel_id2idx)
            df['general_relation'] = df['general_relation'].replace(self.embedding_gen_rel_id2idx)
            
        return dfs 
    
    # function for extracting true facts from the graph and generating negative facts
    def getFacts(self):
        
        # get all true tail entities for each head/relation combo
        self.head_rels = pd.DataFrame(self.embedding_facts.groupby(
            ['head', 'general_relation', 'specific_relation'])['tail'].agg(list))
        self.head_rels['neg_tail'] = self.head_rels.apply(
            lambda row: np.setdiff1d(self.embedding_entities, 
                                     np.array(row['tail'], dtype = 'int32'), 
                                     assume_unique = True),
            axis = 1
        )
        # same but for head entities for each tail/relation combo
        self.tail_rels = pd.DataFrame(self.embedding_facts.groupby(
            ['tail', 'general_relation', 'specific_relation'])['head'].agg(list))
        self.tail_rels['neg_head'] = self.tail_rels.apply(
            lambda row: np.setdiff1d(self.embedding_entities, 
                                     np.array(row['head'], dtype = 'int32'), 
                                     assume_unique = True),
            axis = 1
        )
    
        # get mean number of tail elements per head and mean number of head elements per tail, for negative sampling
        self.mean_tails_per_head = {(M, r): df['tail'].apply(len).mean() for
                                    (M, r), df in self.head_rels.groupby(level = [1, 2])}
        
        self.mean_heads_per_tail = {(M, r): df['head'].apply(len).mean() for
                                    (M, r), df in self.tail_rels.groupby(level = [1, 2])}