import argparse
import pickle
import tensorflow as tf
import os
import optuna
from utils.embedding import CustomEmbeddingModel

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('json', type = str, help = 'Path to graph json file, generated by Neo4j')
    parser.add_argument('-o', '--outpath', type = str, default = 'custom_trials.tsv', help = 'Filepath where HP search output should be saved')
    parser.add_argument('-t', '--test_size', type = float, default = 0.1, help = 'Proportion of facts allocated to test set')
    parser.add_argument('-v', '--validation_size', type = float, default = 0.1, help = 'Proportion of facts allocated to validation set')
    args = parser.parse_args()
    
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
        
    # initialize model
    model = CustomEmbeddingModel(args.json, state = 123)
    
    # train/test/validation split
    model.split(model.embedding_data,
                test_size = args.test_size, 
                validation_size = args.validation_size, 
                state = 123)
        
    train, validation = model.train_set, model.validation_set
    
    # get all positive and negative facts
    negs_per_pos = 5
    pos_facts = train[['head', 'general_relation', 'specific_relation', 'tail']].to_numpy(dtype = 'int32')
    
    # optuna objective function
    def objective(trial):
        
        tf.random.set_seed(123)
        
        hps = {'size_e': trial.suggest_int('size_e', low = 25, high = 100, step = 25),
               'size_r': trial.suggest_int('size_r', low = 25, high = 100, step = 25),
               'margin': trial.suggest_int('margin', low = 0, high = 3, step = 1), 
               'norm': trial.suggest_categorical('norm', [1, 2]),
               'batch_size': trial.suggest_int('batch_size', low = 1000, high = 4000, step = 1000)
               }
        
        lr = trial.suggest_float('lr', low = 0.0001, high = 0.1, log = True)
        
        if hasattr(model, 'optimizer'):
            model.initializeModel(model.embedding_facts, **hps, num_neg_per_pos = negs_per_pos)
            model.optimizer.learning_rate = lr
        else:
            optimizer = tf.keras.optimizers.Adam(learning_rate = lr)
            model.initializeModel(model.embedding_facts, **hps, optimizer = optimizer, num_neg_per_pos = negs_per_pos)
        
        batches = model.getBatches(pos_facts)
        model.train(batches, verbose = False)
        
        amri = model.AMRIEvaluation(validation)
        
        return amri
    
    study = optuna.create_study(direction = 'maximize', 
                                pruner = optuna.pruners.MedianPruner(), 
                                sampler = optuna.samplers.TPESampler(seed = 123))
    
    study.optimize(objective, n_trials = 50)
    
    trials = study.trials_dataframe()
    
    trials.to_csv(args.outpath, sep = '\t', index = False)