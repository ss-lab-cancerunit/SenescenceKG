library(dplyr, quietly = T)
library(stringr, quietly = T)
library(ggplot2, quietly = T)
library(reshape2, quietly = T)

# loop through experiments in the mounted folder
all_data <- list()
metadata <- data.frame()
for(exp_id in dir('data')){
  
  if(grepl('SRP', exp_id)){
    if('counts' %in% dir(paste0('data/', exp_id))){
      all_files <- list.files(paste0('data/', exp_id, '/counts'))
      DE_files <- all_files[grepl('allGene_DESeq2', all_files)]
      DE_filepaths <- paste0('data/', exp_id, '/counts/', DE_files)
      for(file in DE_filepaths){
        res <- read.table(file, header = T)
        res <- dplyr::select(res, GeneID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
        res$GeneID <- str_replace(res$GeneID, '\\..+$', '')
        res <- res[!duplicated(res$GeneID),]
        comparison <- str_extract(file, '(?<=DESeq2_)(.+)(?=\\.csv$)')
        info <- paste(exp_id, comparison, sep = '_')
        colnames(res)[-1] <- paste(colnames(res)[-1], info, sep = '_')
        metadata <- rbind(metadata, c('experiment' = exp_id, 'comparison' = comparison))
        all_data[[info]] <- res
      }
    }
  }
}

colnames(metadata) <- c('experiment', 'comparison')

for(i in 1:length(all_data)){
  
  # get experiment ID and senescence type
  info <- names(all_data)[i]
  
  if(i == 1){
    all_DE_info <- all_data[[info]]
  } else {
    all_DE_info <- dplyr::full_join(all_DE_info, all_data[[info]], 
                                    by = 'GeneID')
  }

}

if(!interactive()){
  write.csv(metadata, file = 'metadata/RNAseq_meta_analysis_metadata_incomplete.csv', row.names = FALSE)
  write.csv(all_DE_info, file = 'data/all_DE_info.csv', row.names = FALSE)  
}
