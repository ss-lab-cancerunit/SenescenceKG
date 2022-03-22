library(metafor)
library(edgeR)
library(ggplot2)
library(stringr)

# performs TMM normalisation and gene filtering
preprocessCounts <- function(counts, samples, genes, design_matrix,
                             min.count = 20, min.total.count = 200,
                             weighting = TRUE){
  
  # initialise DGEList object for edgeR
  DGE <- DGEList(counts = counts, 
                 samples = samples, genes = genes)
  
  # normalise by library size
  normDGE <- calcNormFactors(DGE, method = 'TMM',
                             doWeighting = weighting)
  
  # filter genes with low expression
  genes_to_keep <- filterByExpr(normDGE, 
                                design = design_matrix,
                                min.count = min.count, 
                                min.total.count = min.total.count)
  normDGE$counts <- normDGE$counts[genes_to_keep,]
  normDGE$genes <- normDGE$genes[genes_to_keep,]
  normDGE$filtered_genes <- normDGE$genes[!genes_to_keep,]
  
  return(normDGE)
  
}


# fit a mixed model for a single gene's log fold changes
fitMM <- function(LFCs, weights, FE_form, RE_form, df){
  model <- rma.mv(LFCs, weights, mods = FE_form, random = RE_form, 
                  data = df, test = 't')
  return(model)
}

# conduct differential expression analysis, analysing LFCs from all genes in all studies
lfcMetaAnalysis <- function(genes, LFC_data, SE_data, metadata, FE_form, RE_form){
  
  # initialised results
  results <- list()
  
  # loop through genes, fit mixed model for each
  for(i in 1:nrow(LFC_data)){
    lfc <- unlist(LFC_data[i,])
    SE <- unlist(SE_data[i,])
    gene <- genes[i]
    keep_mask <- which(!is.na(lfc))
    # add LFC into data
    data <- metadata[keep_mask,]
    LFCs <- lfc[keep_mask]
    # get sampling variances from the SEs
    variances <- SE[keep_mask]^2
    # fit mixed model
    if(length(unique(data$senescence)) == length(unique(metadata$senescence))){
      mm_res <- fitMM(LFCs, variances, FE_form, RE_form, data)
      res <- c(mm_res$beta, mm_res$QMp, mm_res$pval) 
      results[[gene]] <- res
    }
    else{
      results[[gene]] <- rep(NA, 7)
    }
  }
  
  results_df <- do.call(rbind, results)
  
  return(results_df)
  
}

# run a permutation test for the gene ranks
permutationTest <- function(gene_rank_matrix, genes_per_study, 
                            sample_size = 100000, state = 123){
  
  set.seed(state)
  
  # sample ranks 
  sampled_ranks <- sapply(genes_per_study, function(n){
    sample(1:n, size = sample_size, replace = T) 
  })
  
  sampled_rrs <- 1 / sampled_ranks
  
  # get gene p values
  gene_pvalues <- apply(gene_rank_matrix, 1, function(generanks){
    
    # filter NA ranks
    valid_rank_mask <- !is.na(generanks)
    valid_ranks <- generanks[valid_rank_mask]
    gene_mrr <- mean(1 / valid_ranks)
    sampled_mrrs <- apply(sampled_rrs[,valid_rank_mask], 1, mean)
    p <- sum(sampled_mrrs <= gene_mrr) / sample_size
    return(p)
    
  })
  
  names(gene_pvalues) <- rownames(gene_rank_matrix)
  
  return(gene_pvalues)
  
}

makeVolcanoPlot <- function(results, 
                            adj_p_column,
                            LFC_column, 
                            title, 
                            LFC_threshold = 1, 
                            p_threshold = 0.05){
  
  results$DE <- ifelse(results[[adj_p_column]] < p_threshold & abs(results[[LFC_column]]) > LFC_threshold, 
                       'Yes', 'No')
  vlc_plot <- ggplot(results, aes(x = !!sym(LFC_column), y = -log10(!!sym(adj_p_column)), 
                                  colour = DE)) +
    geom_point(size = 0.5) + 
    geom_segment(mapping = aes(x = -LFC_threshold, xend = -LFC_threshold,
                               y = 0, yend = Inf), linetype = 'dashed',
                 colour = 'black') +
    geom_segment(mapping = aes(x = LFC_threshold, xend = LFC_threshold,
                               y = 0, yend = Inf), linetype = 'dashed',
                 colour = 'black') +
    geom_segment(mapping = aes(x = -Inf, xend = Inf,
                               y = -log10(p_threshold), yend = -log10(p_threshold)), linetype = 'dashed',
                 colour = 'black') +
    theme_bw() +
    labs(x = '\nEstimated log2 fold change', y = '-log10(adjused p-value)\n', 
         title = title, colour = 'Differentially-expressed') +
    theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(override.aes = list(size = 10)))
  
  return(vlc_plot)
  
}

# function for looking up LFCs
LFClookup <- function(genes, 
                      siggenelist,
                      res){
  
  LFCs <- lapply(genes, function(gene){
    
    # get the senescence types identified as being significantly associated with each gene
    senescence_type_mask <- sapply(siggenelist, function(genes) gene %in% genes)
    if(any(senescence_type_mask)){
      senescence_types <- str_flatten(names(siggenelist)[senescence_type_mask], collapse = '/')
    } else {
      senescence_types <- 'None'
    }
    
    # lookup LFCs from p-value tables (assumes rownames are genes)
    cols <- c('LFC_Oncogene', 'LFC_Replicative', 'LFC_DNAdamage',
              'adj_p_Oncogene', 'adj_p_Replicative', 'adj_p_DNAdamage')
    gene_data <- res[gene, cols]
    gene_data['sig_types'] <- senescence_types
    
    return(gene_data)
    
  })
  
  LFC_df <- do.call(rbind, LFCs)
  rownames(LFC_df) <- genes
  
  return(LFC_df)
  
}
