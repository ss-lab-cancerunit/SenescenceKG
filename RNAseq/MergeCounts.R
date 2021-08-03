library(dplyr, quietly = T)
library(stringr, quietly = T)
library(ggplot2, quietly = T)
library(reshape2, quietly = T)
library(edgeR, quietly = T)
library(org.Hs.eg.db, quietly = T)
library(tibble, quietly = T)

source('utils.R')

# loop through experiments in the mounted folder
all_counts <- list()
for(exp_id in dir('data')){
  
  if(grepl('SRP', exp_id)){
    if('counts' %in% dir(paste0('data/', exp_id))){
      # load rdata file with counts
      filepath <- paste0('data/', exp_id, '/counts/', exp_id, '_experiment.json_counts.RData')
      load(file = filepath)
      # add experiment id to column names
      colnames(counts$counts) <- paste(colnames(counts$counts), exp_id, sep = '_')
      # save gene IDs, convert to tibble, add column for gene
      gene_ids <- stringr::str_replace(rownames(counts$counts), '\\..+$', '')
      counts$counts <- as_tibble(counts$counts)
      counts$counts$gene <- gene_ids
      counts$counts <- relocate(counts$counts, gene)
      counts$counts <- dplyr::group_by(counts$counts, gene) %>%
        dplyr::summarise(across(.fns = sum))
      # save to list
      all_counts[[exp_id]] <- counts
    }
  }
}

# merge counts into single dataframe
for(i in 1:length(all_counts)){
  
  # get experiment ID and senescence type
  exp_id <- names(all_counts)[i]
  
  if(i == 1){
    count_matrix <- all_counts[[exp_id]]$counts
  } else {
    count_matrix <- dplyr::full_join(count_matrix, all_counts[[exp_id]]$counts, 
                                     by = 'gene', all = TRUE)
  }
  
}

count_matrix[is.na(count_matrix)] <- 0
count_columns <- colnames( dplyr::select(count_matrix, -gene) )
# make an initial metadata file
sample_metadata <- data.frame('variable' = count_columns,
                              'experiment' = str_extract(count_columns, 'SRP.+$'), 
                              'senescent' = ifelse(grepl('Control', count_columns), 
                                                   'Control', 'Senescent'))
# if run as a script, save files
if(!interactive()){
  write.csv(sample_metadata, file = 'metadata/RNAseq_sample_metadata_incomplete.csv', row.names = FALSE)
  write.csv(count_matrix, file = 'data/combined_counts.csv', row.names = FALSE) 
}

