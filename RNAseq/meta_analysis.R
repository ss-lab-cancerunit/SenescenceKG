library(metafor)
library(dplyr)
library(tibble)

source('utils.R')
source('merge_DE.R')

# extract log fold changes and standard errors
all_lfcs <- select_if(all_DE_info, grepl('log2FoldChange', colnames(all_DE_info)))
keep <- apply(all_lfcs, 1, function(row) sum(is.na(row)) < 10)
keep_lfcs <- all_lfcs[keep,]

all_SEs <- select_if(all_DE_info, grepl('lfcSE', colnames(all_DE_info)))
keep_SEs <- all_SEs[keep,]

keep_genes <- all_DE_info[keep, 'GeneID']

# import metadata
metadata <- read.csv('metadata/RNAseq_meta_analysis_metadata.csv')
metadata$senescence <- factor(metadata$senescence, levels = c('DNA-damage', 'Oncogene', 'Replicative'))

# set model formula
FE_form <- ~ senescence - 1 
RE_form <- ~ (1 | experiment)

writeLines('Running mixed model meta-analysis...')

# run meta analysis
MM_results <- lfcMetaAnalysis(keep_genes,
                              keep_lfcs, keep_SEs, 
                              metadata, FE_form, RE_form)
MM_results <- as.data.frame(MM_results)
colnames(MM_results) <- c('LFC_DNAdamage', 
                          'LFC_Oncogene',
                          'LFC_Replicative',
                          'F_test', 
                          'p_DNAdamage', 
                          'p_Oncogene',
                          'p_Replicative')
MM_results <- rownames_to_column(MM_results, 'ensembl')

# adjust p values
adjusted_p_values <- p.adjust(c(MM_results$p_DNAdamage,
                                MM_results$p_Oncogene,
                                MM_results$p_Replicative),
                              method = 'BH')
adj_p_matrix <- matrix(adjusted_p_values, ncol = 3)
MM_results[,c('adj_p_DNAdamage', 'adj_p_Oncogene', 'adj_p_Replicative')] <- adj_p_matrix

write.csv(MM_results, 'data/meta_analysis_results.csv', row.names = F)


