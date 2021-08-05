import numpy as np
import tensorflow as tf

class TransR:
    
    def __init__(self, 
                 size_e: int, 
                 size_r: int, 
                 num_e: int, 
                 num_r: int,
                 num_m: int,
                 norm: int,
                 margin: int, 
                 initializer: tf.keras.initializers,
                 optimizer: tf.keras.optimizers):
        
        assert norm in [1, 2], 'Norm should be either 1 or 2'
        
        self.size_e = size_e
        self.size_r = size_r
        self.num_e = num_e
        self.num_r = num_r
        self.num_m = num_m
        self.norm = norm
        self.margin = margin
        self.optimizer = optimizer
        
        # initialise entity/relation embeddings with random values
        self.ent_embeddings = tf.Variable(initializer(shape = [self.num_e, self.size_e], dtype = tf.float32),
                                          name = 'ent_embedding', dtype = tf.float32)
                    
        self.rel_embeddings = tf.Variable(initializer(shape = [self.num_r, self.size_r], dtype = tf.float32),
                                          name = 'rel_embedding', dtype = tf.float32)
        
            
        # l2 normalisation
        self.ent_embeddings.assign(tf.math.l2_normalize(self.ent_embeddings, axis = 1))
        self.rel_embeddings.assign(tf.math.l2_normalize(self.rel_embeddings, axis = 1))
        
        # initialise relation matrices 
        self.rel_matrices = tf.Variable(initializer(shape = [self.num_m, self.size_r, self.size_e], dtype = tf.float32), 
                                        name = 'rel_matrices', dtype = tf.float32)
            
    
    @property
    def params(self):
        return [self.ent_embeddings, self.rel_embeddings, self.rel_matrices]
    
    # fact scoring function, takes 4 vectors as input: head (h), taiL (t), general relation (M), and specific relation (r)
    def score(self,
              h: np.ndarray, 
              t: np.ndarray, 
              r: np.ndarray, 
              M: np.ndarray):
        
        # lookup embeddings for entities and relations
        h_e = tf.nn.embedding_lookup(self.ent_embeddings, h)
        t_e = tf.nn.embedding_lookup(self.ent_embeddings, t)
        r_e = tf.nn.embedding_lookup(self.rel_embeddings, r)
        M_e = tf.nn.embedding_lookup(self.rel_matrices, M)
        
        # project entity embeddings with M matrix
        h_proj = tf.math.l2_normalize(tf.reshape(tf.matmul(M_e, h_e[:, :, None]), (tf.shape(h_e)[0], tf.shape(M_e)[1])), axis = 1)
        t_proj = tf.math.l2_normalize(tf.reshape(tf.matmul(M_e, t_e[:, :, None]), (tf.shape(t_e)[0], tf.shape(M_e)[1])), axis = 1)
        
        # calculate score
        if self.norm == 2:
            score = tf.math.reduce_sum((h_proj + r_e - t_proj)**2, axis = 1)
        else:
            score = tf.math.reduce_sum(abs(h_proj + r_e - t_proj), axis = 1)
            
        return score
            
    # calculate the margin loss for a batch of facts
    def marginLoss(self, 
                   pos_h: np.ndarray, 
                   pos_t: np.ndarray, 
                   neg_h: np.ndarray, 
                   neg_t: np.ndarray, 
                   pos_r: np.ndarray,
                   neg_r: np.ndarray,
                   pos_M: np.ndarray,
                   neg_M: np.ndarray):
        
        # get scores for positive and negative entities
        pos_scores = self.score(pos_h, pos_t, pos_r, pos_M)
        neg_scores = self.score(neg_h, neg_t, neg_r, neg_M)
            
        # calculate total loss
        loss = tf.math.reduce_sum(tf.maximum(pos_scores[:,None] - neg_scores[None,:] + self.margin, 0))
            
        return loss
  
    def trainStep(self, 
                  pos_examples: np.ndarray, 
                  neg_examples: np.ndarray):
        
        # initialise gradient tracker, calculate loss
        with tf.GradientTape() as tape:
            
            # normalise embeddings
            self.ent_embeddings.assign(tf.math.l2_normalize(self.ent_embeddings, axis = 1), name = 'ent_embedding')
            self.rel_embeddings.assign(tf.math.l2_normalize(self.rel_embeddings, axis = 1), name = 'rel_embedding')
            
            # split into head, general relation (M), specific relation (r), and tail vectors
            pos_h, M_pos, r_pos, pos_t = pos_examples[:,0], pos_examples[:,1], pos_examples[:,2], pos_examples[:,3]
            neg_h, M_neg, r_neg, neg_t = neg_examples[:,0], neg_examples[:,1], neg_examples[:,2], neg_examples[:,3]
             
            # calculate loss
            loss = self.marginLoss(pos_h, pos_t, neg_h, neg_t, r_pos, r_neg, M_pos, M_neg)
                
        # apply gradients to update the weights
        grads = tape.gradient(loss, self.params)
        self.optimizer.apply_gradients(zip(grads, self.params))
        
        return loss