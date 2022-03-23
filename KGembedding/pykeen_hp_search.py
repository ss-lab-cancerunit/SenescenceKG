import pandas as pd
import argparse
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory
from pykeen.hpo import hpo_pipeline
from torch.nn.init import xavier_normal_
import numpy as np
import pickle
import warnings
warnings.filterwarnings('ignore')

from utils.preprocessing import GraphParser

if __name__ == '__main__':

    ## Arguments for the script are the graph json file and output directory
    parser = argparse.ArgumentParser()
    parser.add_argument('json', type = str, help = 'Path to graph json file, generated by Neo4j')
    parser.add_argument('-m', '--model', type = str, help = 'Which PyKeen embedding model to use (as string)')
    parser.add_argument('-o', '--outpath', type = str, default = 'trials.tsv', help = 'Filepath where hyperparameter search results should be saved')
    parser.add_argument('-t', '--test_size', type = float, default = 0.1, help = 'Proportion of facts allocated to test set')
    parser.add_argument('-v', '--validation_size', type = float, default = 0.1, help = 'Proportion of facts allocated to validation set')
    args = parser.parse_args()

    # load the file into the graph parser class
    graphdata = GraphParser(args.json)
    
    # filter relations that occur at low frequency
    graphdata.processGraph()
        
    # train/test/validation split
    facts = graphdata.embedding_data.drop_duplicates(['head_lab', 'general_relation_lab', 'tail_lab'])
    graphdata.split(facts, 
                    test_size = args.test_size, 
                    validation_size = args.validation_size, 
                    stratify_column = 'general_relation_lab',
                    state = 123)
    train = graphdata.train_set[['head_lab', 'general_relation_lab', 'tail_lab']].to_numpy(dtype = 'str')
    test = graphdata.test_set[['head_lab', 'general_relation_lab', 'tail_lab']].to_numpy(dtype = 'str')
    validation = graphdata.validation_set[['head_lab', 'general_relation_lab', 'tail_lab']].to_numpy(dtype = 'str')
    
    # make pykeen triples objects
    train_triples = TriplesFactory.from_labeled_triples(train, 
                                                        create_inverse_triples = False,
                                                        entity_to_id = graphdata.node_lab2id, 
                                                        relation_to_id = graphdata.embedding_gen_rel_lab2idx)
    validation_triples = TriplesFactory.from_labeled_triples(validation, 
                                                             create_inverse_triples = False,
                                                             entity_to_id = graphdata.node_lab2id, 
                                                             relation_to_id = graphdata.embedding_gen_rel_lab2idx)
    test_triples = TriplesFactory.from_labeled_triples(test,
                                                       create_inverse_triples = False,
                                                       entity_to_id = graphdata.node_lab2id, 
                                                       relation_to_id = graphdata.embedding_gen_rel_lab2idx)
    
    # add additional parameter to search through for transR models
    if args.model.lower() == 'transr':
        model_hp_ranges = {'embedding_dim': dict(type = int, low = 25, high = 100, step = 25),
                           'relation_dim': dict(type = int, low = 25, high = 100, step = 25)}
    else:
        model_hp_ranges = {'embedding_dim': dict(type = int, low = 25, high = 100, step = 25)}
    
    # build the model
    results = hpo_pipeline( 
        n_trials = 50,
        training = train_triples,
        testing = test_triples,
        validation = validation_triples,
        model = args.model,
        training_kwargs = {'num_epochs': 50},
        training_kwargs_ranges = {'batch_size': dict(type = int, low = 1000, high = 4000, step = 1000)},
        model_kwargs = {'entity_initializer': xavier_normal_, 'relation_initializer': xavier_normal_},
        model_kwargs_ranges = model_hp_ranges,
        negative_sampler_kwargs = {'num_negs_per_pos': 5},
        loss = 'MarginRankingLoss',
        optimizer = 'Adam',
        optimizer_kwargs_ranges = {'lr': dict(type = float, low = 0.0001, high = 0.1, log = True)},
        negative_sampler = 'Bernoulli',
        training_loop = 'sLCWA',
        sampler_kwargs = {'seed': 123},
        stopper = 'early'
    )
    
    # save
    trials = results.study.trials_dataframe()
    trials.to_csv(args.outpath, sep = '\t', index = False)