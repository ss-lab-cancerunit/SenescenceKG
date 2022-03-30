import pandas as pd
import numpy as np 
import pickle
from tensorflow import keras
import torch
import optuna
import tensorflow as tf
import xml.etree.ElementTree as ET
import logging as lg
import sys
from sklearn.model_selection import StratifiedKFold, train_test_split

logger = lg.getLogger(__name__)

formatter = lg.Formatter('%(asctime)s %(levelname)-8s: %(message)s')
streamhandler = lg.StreamHandler(sys.stdout)
streamhandler.setLevel(lg.INFO)
streamhandler.setFormatter(formatter)
logger.addHandler(streamhandler)
logger.setLevel(lg.DEBUG)
logger.propagate = False

# make class for optuna objective

class Objective:
    
    def __init__(self, X_train, y_train):
        
        self.X_train = X_train
        self.y_train = y_train
        y_split = y_train.argmax(axis = 1) if len(self.y_train.shape) > 1 else y_train
        self.X_train_folds = list(StratifiedKFold(n_splits = 5).split(X_train, y_split))
        self.optimizer = keras.optimizers.Adam()
        
    def __call__(self, trial):
        
        np.random.seed(123)
        tf.random.set_seed(123)
        
        lr = trial.suggest_float('lr', low = 0.0001, high = 0.1, log = True)
        hidden_layer_size = trial.suggest_int('hidden_layer_size', low = 20, high = 100, step = 10)
        activation_func = trial.suggest_categorical('activation', ['relu', 'sigmoid', 'tanh'])
        l2_regularization_weight = trial.suggest_float('l2_regularization', low = 0, high = 0.1, step = 0.01)
        hidden_dropout_rate = trial.suggest_float('hidden_dropout_rate', low = 0, high = 0.5, step = 0.1)
        num_epochs = trial.suggest_int('num_epochs', low = 10, high = 50, step = 5)
        
        metric_values = []
        
        final_layer_activation = 'softmax' if len(self.y_train.shape) > 1 else 'sigmoid'
        final_layer_size = self.y_train.shape[1] if len(self.y_train.shape) > 1 else 1
        loss_func = 'categorical_crossentropy' if len(self.y_train.shape) > 1 else 'binary_crossentropy'
        metric = 'accuracy' if len(self.y_train.shape) > 1 else keras.metrics.AUC()
        
        # get mean accuracy across all 5 folds of training data
        for train, validation in self.X_train_folds:
            
            model = keras.Sequential([keras.layers.InputLayer(input_shape = (self.X_train.shape[1],)), 
                                      keras.layers.Dense(hidden_layer_size, activation = activation_func, 
                                                         kernel_regularizer = keras.regularizers.l2(l2_regularization_weight),
                                                         kernel_initializer = 'glorot_normal'),
                                      keras.layers.Dropout(hidden_dropout_rate),
                                      keras.layers.Dense(final_layer_size, 
                                                         activation = final_layer_activation, 
                                                         kernel_regularizer = keras.regularizers.l2(l2_regularization_weight),
                                                         kernel_initializer = 'glorot_normal')])
        
            self.optimizer.lr = lr
            
            model.compile(loss = loss_func, 
                          optimizer = self.optimizer,
                          metrics = [metric])
            
            X = self.X_train[train]
            y = self.y_train[train]
            model.fit(X, y, epochs = num_epochs, batch_size = len(X))
            error, value = model.evaluate(self.X_train[validation], self.y_train[validation])
            metric_values.append(value)
                        
        mean_value = float(np.mean(metric_values))
        
        return mean_value


if __name__ == '__main__':
    
    # import models, get embeddings

    transE_model = pickle.load(open('data/embeddings/transE_model.p', 'rb'))
    transE_embeddings = transE_model['model'].entity_representations[0].cpu()
    transR_model = pickle.load(open('data/embeddings/transR_model.p', 'rb'))
    transR_embeddings = transR_model['model'].entity_representations[0].cpu()
    convE_model = pickle.load(open('data/embeddings/convE_model.p', 'rb'))
    convE_embeddings = convE_model['model'].entity_representations[0].cpu()
    complEx_model = pickle.load(open('data/embeddings/complEx_model.p', 'rb'))
    complEx_embeddings = complEx_model['model'].entity_representations[0].cpu()
    custom_model = pickle.load(open('data/embeddings/custom_model.p', 'rb'))
    custom_embeddings = custom_model['ent_embeddings']

    logger.info(f'Embedding models imported.')
    
    # import train/test sets used in the embedding models
    
    pykeen_train = pd.read_csv('data/embeddings/pykeen_train.tsv', sep = '\t')
    custom_train = pd.read_csv('data/embeddings/custom_train.tsv', sep = '\t')
    
    # make dictionaries mapping entity index to entity id 
    
    custom_all_entities = custom_train['head'].append(custom_train['tail'])
    custom_all_entity_ids = custom_train['head_id'].append(custom_train['tail_id'])
    custom_entity_id2idx = {ent_id: ent_idx for ent_idx, ent_id in zip(custom_all_entities, custom_all_entity_ids)}
    
    pykeen_all_entities = pykeen_train['head'].append(pykeen_train['tail'])
    pykeen_all_entity_ids = pykeen_train['head_id'].append(pykeen_train['tail_id'])
    pykeen_entity_id2idx = {ent_id: ent_idx for ent_idx, ent_id in zip(pykeen_all_entities, pykeen_all_entity_ids)}
    
    # extract gene IDs from the node dictionary
    
    gene_entities = pd.Series([node_id for node_id, entry in custom_model['entities'].items() if entry['labels'] == ['Gene']])
    
    # get gene entities from pykeen and custom embeddings
    
    pykeen_gene_ids = set([ent for ent in pykeen_all_entity_ids
                           if ent in gene_entities.values])
    custom_gene_ids = set([ent for ent in custom_all_entity_ids
                           if ent in gene_entities.values])
    
    # subset gene entities that appear in both pykeen and custom model
    gene_entities = pd.Series(list(pykeen_gene_ids & custom_gene_ids))
    
    # get genes with trained embeddings
    valid_LFC_mask = np.array([custom_model['entities'][gene]['properties']['replicative_LFC'] != 'NA' and
                               custom_model['entities'][gene]['properties']['oncogene_LFC'] != 'NA' and 
                               custom_model['entities'][gene]['properties']['DNAdamage_LFC'] != 'NA' for gene in gene_entities])
    valid_genes = gene_entities[valid_LFC_mask]
    valid_gene_ensembl_ids = [custom_model['entities'][gene]['properties']['ensembl'] for gene in valid_genes]
    
    pykeen_gene_idx = [pykeen_entity_id2idx[gene] for gene in valid_genes]
    custom_gene_idx = [custom_entity_id2idx[gene] for gene in valid_genes]
    
    # extract gene embeddings
    
    gene_id_tensor = torch.tensor(pykeen_gene_idx, dtype = torch.long)
    
    embeddings = {'transE': transE_embeddings(gene_id_tensor).detach().numpy(), 
                  'transR': transR_embeddings(gene_id_tensor).detach().numpy(), 
                  'convE': convE_embeddings(gene_id_tensor).detach().numpy(), 
                  'complEx': complEx_embeddings(gene_id_tensor).detach().numpy(), 
                  'custom': tf.gather(custom_embeddings, custom_gene_idx).numpy()
                 }
    
    ## get target values for each gene (senescence subtypes and pathways)
    
    # make array of classifications for all 3 subtypes of senescence
    replicative_senescence_target = np.zeros(shape = (len(valid_genes), 3))
    oncogene_senescence_target = np.zeros(shape = (len(valid_genes), 3))
    DNAdamage_senescence_target = np.zeros(shape = (len(valid_genes), 3))
    
    # loop through genes and add their labels to the target arrays
    for i, gene in enumerate(valid_genes):
        replicative_LFC = custom_model['entities'][gene]['properties']['replicative_LFC']
        if replicative_LFC > 1:
            replicative_senescence_target[i, 2] = 1
        elif replicative_LFC < -1:
            replicative_senescence_target[i, 0] = 1
        else:
            replicative_senescence_target[i, 1] = 1
        
        oncogene_LFC = custom_model['entities'][gene]['properties']['oncogene_LFC']
        if oncogene_LFC > 1:
            oncogene_senescence_target[i, 2] = 1
        elif oncogene_LFC < -1:
            oncogene_senescence_target[i, 0] = 1
        else:
            oncogene_senescence_target[i, 1] = 1
            
        DNAdamage_LFC = custom_model['entities'][gene]['properties']['DNAdamage_LFC']
        if DNAdamage_LFC > 1:
            DNAdamage_senescence_target[i, 2] = 1
        elif DNAdamage_LFC < -1:
            DNAdamage_senescence_target[i, 0] = 1
        else:
            DNAdamage_senescence_target[i, 1] = 1
            
    # import gene sets for the three pathways to be predicted
    geneset_filepaths = ['data/genesets/SASP_genes.xml', 
                         'data/genesets/Cell_cycle_checkpoint_genes.xml',
                         'data/genesets/TP53_regulation_genes.xml']
    geneset_names = ['SASP', 'Cell cycle checkpoint', 'TP53 regulation']
    
    # put all targets into one dictionary
    targets = {'replicative': replicative_senescence_target,
               'oncogene': oncogene_senescence_target,
               'DNA damage': DNAdamage_senescence_target,
              }
    
    # parse MSigDB's xml files and check if each gene's ensembl id is present in each set
    for filepath, name in zip(geneset_filepaths, geneset_names):
        
        tree = ET.parse(filepath)
        root = tree.getroot()
        pathway_genes = root.find('GENESET').get('MEMBERS').split(',')
        targets[name] = np.array([int(gene in pathway_genes) for gene in valid_gene_ensembl_ids])
    
    # get train and test indices
    train, test = train_test_split(np.arange(len(valid_genes)), test_size = 0.2, random_state = 123)
    
    # run optuna optimization for each target with each model's embeddings
    results = {'train': train, 'test': test, 'gene_ids': valid_genes, 
               'custom_model_id2idx': custom_entity_id2idx,
               'pykeen_model_id2idx': pykeen_entity_id2idx,
               'embeddings': embeddings, 'targets': targets}
    
    for i, target in enumerate(targets, start=1):
        y = targets[target]
        results[target] = {}
        for model in embeddings:
            X = embeddings[model]
            objective = Objective(X[train], y[train])
            study = optuna.create_study(direction = 'maximize', 
                                        pruner = optuna.pruners.MedianPruner(), 
                                        sampler = optuna.samplers.TPESampler(seed = 123))
            study.optimize(objective, n_trials = 50)
            trials = study.trials_dataframe()
            results[target][model] = trials
        logger.info(f'Hyperparameter search for target {i} completed')
            
    pickle.dump(results, open('data/classifier/classifier_hp_results.p', 'wb'))