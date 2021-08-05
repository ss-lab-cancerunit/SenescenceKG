import numpy as np
import pandas as pd
import json
import tensorflow as tf
import optuna
import random
from numba import jit
from scipy.stats import rankdata

from utils.transR import TransR
from utils.preprocessing import GraphParser

class CustomEmbeddingModel(TransR, GraphParser):
    
    def __init__(self, 
                 json_filepath: str,
                 state: int = 123,
                 *args, **kwargs):
        
        GraphParser.__init__(self, json_filepath)
        
        self.processGraph(*args, **kwargs)
        
        self.state = state
        
        print('Graph data parsed and cleaned.\n')
         
    def initializeModel(self,
                        facts: pd.DataFrame,
                        size_e: int = 50,
                        size_r: int = 25,
                        margin: float = 1,
                        norm: int = 2, 
                        batch_size: int = 1000,
                        num_neg_per_pos: int = 5,
                        n_epochs: int = 50,
                        initializer_class: tf.keras.initializers = None,
                        optimizer: tf.keras.optimizers = None
                       ):
        
        if not hasattr(self, 'initializer') and initializer_class is None:
            self.initializer = tf.keras.initializers.glorot_normal(seed = self.state)
        elif initializer_class is not None:
            self.initializer = initializer_class(seed = self.state)
            
        if not hasattr(self, 'optimizer') and optimizer is None:
            self.optimizer = tf.keras.optimizers.Adam()
        elif optimizer is not None:
            self.optimizer = optimizer
            
        num_entities = len(pd.unique(facts['head'].append(facts['tail'])))
        num_specific_relations = len(pd.unique(facts['specific_relation']))
        num_general_relations = len(pd.unique(facts['general_relation']))
    
        TransR.__init__(self, 
                        size_e = size_e,
                        size_r = size_r,
                        num_e = num_entities,
                        num_r = num_specific_relations,
                        num_m = num_general_relations,
                        margin = margin,
                        norm = norm,
                        initializer = self.initializer,
                        optimizer = self.optimizer)
        
        self.batch_size = batch_size
        self.num_neg_per_pos = num_neg_per_pos
        self.n_epochs = n_epochs
        
    # negative sampling of entities
    @tf.function
    def corruptTriplets(self,
                        facts: tf.Tensor,
                        p_head: np.ndarray,
                        n: int):
        
        # get head/tail indices to replace
        replace_head = tf.random.uniform(shape = [len(p_head) * n], minval = 0, 
                                         maxval = 1, dtype = tf.float32) < tf.repeat(p_head, n)
        replace_tail = ~replace_head
        replace_head = tf.reshape(tf.where(replace_head), [-1])
        replace_tail = tf.reshape(tf.where(replace_tail), [-1])
         
        # sample negative entities
        sampled_negs = tf.random.uniform(shape = [len(facts) * n], minval = 0, 
                                         maxval = self.num_entities - 1, 
                                         dtype = tf.int32)
        
        # make an array of negative facts
        facts_expanded = tf.repeat(facts, n, axis = 0)
        
        # identify where sampled entities match entities in true facts
        heads_match = tf.gather(sampled_negs, replace_head) == tf.gather(facts_expanded[:, 0], replace_head)
        tails_match = tf.gather(sampled_negs, replace_tail) == tf.gather(facts_expanded[:, 3], replace_tail)
        
        # if the sampled entities match true facts, add 1
        new_heads = tf.where(heads_match,
                             tf.gather(sampled_negs, replace_head) + tf.cast(heads_match, dtype = tf.int32),
                             tf.gather(sampled_negs, replace_head))
        new_tails = tf.where(tails_match,
                             tf.gather(sampled_negs, replace_tail) + tf.cast(tails_match, dtype = tf.int32),
                             tf.gather(sampled_negs, replace_tail))
        
        # make corrupted head and tail fact tensors
        corrupted_head_facts = tf.transpose([new_heads, 
                                             tf.gather(facts_expanded[:, 1], replace_head),
                                             tf.gather(facts_expanded[:, 2], replace_head),
                                             tf.gather(facts_expanded[:, 3], replace_head)]
                                            )
        corrupted_tail_facts = tf.transpose([tf.gather(facts_expanded[:, 0], replace_tail), 
                                             tf.gather(facts_expanded[:, 1], replace_tail),
                                             tf.gather(facts_expanded[:, 2], replace_tail),
                                             new_tails]
                                            )
        
        # put into a single tensor and sort to match the order of the input facts
        all_neg_facts = tf.concat([corrupted_head_facts, corrupted_tail_facts], 0)
        sorted_index = tf.argsort(tf.concat([replace_head, replace_tail], 0))
        sorted_neg_facts = tf.gather(all_neg_facts, sorted_index) 

        return sorted_neg_facts
        
    # function for corrupting facts via bernoulli sampling
    def bernoulliSample(self,
                        facts: np.ndarray, 
                        n: int) -> tf.Tensor:
        
        # get the mean number of tail entities associated with each head for general relation M/specific relation r
        mean_tails_per_head = np.empty(shape = len(facts), dtype = 'float32')
        # get the mean number of head entities associated with each tail for general relation M/specific relation r
        mean_heads_per_tail = np.empty(shape = len(facts), dtype = 'float32')
        for i, rels in enumerate(zip(facts[:,1], facts[:,2])):
            mean_tails_per_head[i] = self.mean_tails_per_head[rels]
            mean_heads_per_tail[i] = self.mean_heads_per_tail[rels]

        # get probability of corrupting head
        p_head = mean_heads_per_tail / (mean_heads_per_tail + mean_tails_per_head)
        
        # get negative facts
        neg_facts = self.corruptTriplets(facts, p_head, n)
        
        return neg_facts
    
    # split a provided data set into batches
    def getBatches(self, 
                   facts: np.ndarray) -> list:
        
        # split data into batches
        batches = []
        for i in range(0, facts.shape[0], self.batch_size):
            start = i
            end = min(start + self.batch_size, facts.shape[0])
            batch = facts[start:end]
            batches.append(batch)
                
        return batches
    
    # train a model on a list of batches
    def train(self, batches: list, verbose: bool  = True, return_results: bool = False):
        
        random.seed(self.state)
     
        for e in range(self.n_epochs):
            
            # initialise total loss to 0
            loss = 0
            
            # shuffle batches
            random.shuffle(batches)
            
            # loop over batches, add to loss
            for pos_facts in batches:
                
                neg_facts = self.bernoulliSample(pos_facts, n = self.num_neg_per_pos)
                loss += self.trainStep(pos_facts, neg_facts) 
                
            if verbose == True:
                print('Loss: {}, epoch {}'.format(loss.numpy(), e+1))
            
        if return_results == True:
            
            results = {'ent_embeddings': self.ent_embeddings,
                       'rel_embeddings': self.rel_embeddings,
                       'rel_matrices': self.rel_matrices,
                       'general_relations': self.gen_rel_id2lab,
                       'specific_relations': self.spec_rel_id2lab,
                       'entities': self.node_entries
                      }
            
            return results
    
    # given a relation h, M, r, t and a set of negative head/tail entities, generate an evaluation set
    @staticmethod
    @jit(nopython = True)
    def makeEvaluationSet(h, M, r, t, neg_h, neg_t) -> np.ndarray:

        facts = np.empty(shape = (1 + len(neg_t) + len(neg_h), 4), dtype = 'int32')
        facts[:(len(neg_t)+1), 0] = h
        facts[(len(neg_t)+1):, 0] = neg_h
        facts[0, 3] = t
        facts[1:(len(neg_t)+1), 3] = neg_t
        facts[(len(neg_t)+1):, 3] = t
        facts[:,1] = M
        facts[:,2] = r
        
        return facts
    
    # calculate an AMRI score for a validation set
    def AMRIEvaluation(self, 
                       validation_set: pd.DataFrame,
                       batch_size: int = 5) -> float:
    
        amri_numerator = 0
        amri_denominator = 0
        
        validation_facts = validation_set[['head', 'general_relation', 'specific_relation', 'tail']].to_numpy(dtype = 'int32')
        
        # loop over batches in the validation set
        for i in range(0, len(validation_facts), batch_size):
            
            max_idx = min(i + batch_size, len(validation_facts))
            
            # get the total size of the evaluation sets for the entire batch of facts
            total_size = sum([len(self.tail_rels.loc[tuple(validation_facts[j, [3, 1, 2]]), 'neg_head']) +
                              len(self.head_rels.loc[tuple(validation_facts[j, [0, 1, 2]]), 'neg_tail']) for
                              j in range(i, max_idx)]) + batch_size
            
            batch_eval_set = np.empty(shape = (total_size, 4), dtype = 'int32')
            
            # fill in the evaluation set array
            start = 0
            fact_evaluation_indices = {}
            for j in range(i, max_idx):
                h, M, r, t = validation_facts[j]
                neg_h = np.array(self.tail_rels.loc[(t, M, r), 'neg_head'])
                neg_t = np.array(self.head_rels.loc[(h, M, r), 'neg_tail'])
                end = start + len(neg_h) + len(neg_t) + 1
                fact_evaluation_indices[j] = start, end
                batch_eval_set[start:end] = self.makeEvaluationSet(h, M, r, t, neg_h, neg_t)
                start = end

            # calculate scores for all facts
            scores = self.score(h = batch_eval_set[:,0],
                                M = batch_eval_set[:,1],
                                r = batch_eval_set[:,2],  
                                t = batch_eval_set[:,3]).numpy()
            
            # get each score's rank among facts in its evaluation set
            for j in range(i, max_idx):
                start, end = fact_evaluation_indices[j]
                ranks = rankdata(scores[start:end], method = 'average')
                amri_numerator += (ranks[0] - 1)
                amri_denominator += len(scores[start:end]) 

        # calculate AMRI
        AMRI = 1 - (2 * amri_numerator / amri_denominator)

        return AMRI
    
    # get head predictions for a general relation, specific relation, and tail
    def getHeadPredictions(self, tail, general_relation, specific_relation):
        
        # get set of possible head entities (not included in graph)
        if (tail, general_relation, specific_relation) in self.tail_rels.index:
            possible_heads = self.tail_rels.loc[(tail, general_relation, specific_relation), 'neg_head']
        else:
            possible_heads = self.embedding_entities
        
        # make data frame of possible facts
        candidate_set = pd.DataFrame({'head': possible_heads,
                                      'general_relation': np.repeat(general_relation, len(possible_heads)),
                                      'specific_relation': np.repeat(specific_relation, len(possible_heads)),
                                      'tail': np.repeat(tail, len(possible_heads))}
                                    )
        
        # score each fact
        candidate_scores = self.score(h = candidate_set['head'].to_numpy(),
                                      M = candidate_set['general_relation'].to_numpy(),
                                      r = candidate_set['specific_relation'].to_numpy(),
                                      t = candidate_set['tail'].to_numpy()).numpy()
        
        # sort data frame by score
        candidate_set['score'] = -candidate_scores
        candidate_set = candidate_set.sort_values('score', ascending = False)
        
        return candidate_set
    
    # get tail predictions for a general relation, specific relation, and head
    def getTailPredictions(self, general_relation, specific_relation, head):
        
        # get set of possible head entities (not included in graph)
        if (head, general_relation, specific_relation) in self.head_rels.index: 
            possible_tails = self.head_rels.loc[(head, general_relation, specific_relation), 'neg_tail']
        else:
            possible_tails = self.embedding_entities
        
        # make data frame of possible facts
        candidate_set = pd.DataFrame({'head': np.repeat(head, len(possible_tails)),
                                      'general_relation': np.repeat(general_relation, len(possible_heads)),
                                      'specific_relation': np.repeat(specific_relation, len(possible_heads)),
                                      'tail': possible_tails}
                                    )
        
        # score each fact
        candidate_scores = self.score(h = candidate_set['head'].to_numpy(),
                                      M = candidate_set['general_relation'].to_numpy(),
                                      r = candidate_set['specific_relation'].to_numpy(),
                                      t = candidate_set['tail'].to_numpy()).numpy()
        
        # sort data frame by score
        candidate_set['score'] = -candidate_scores
        candidate_set = candidate_set.sort_values('score', ascending = False)
        
        return candidate_set