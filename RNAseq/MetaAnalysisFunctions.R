library(metafor)

# performs TMM normalisation and gene filtering
preprocessCounts <- function(counts, samples, genes, design_matrix,
                             min.count = 20, min.total.count = 200){
  
  # initialise DGEList object for edgeR
  DGE <- DGEList(counts = counts, 
                 samples = metadata, genes = genes)
  
  # normalise by library size
  normDGE <- calcNormFactors(DGE, method = 'TMM',
                             doWeighting = FALSE)
  
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
    p <- sum(gene_mrr <= sampled_mrrs) / sample_size
    return(p)
    
  })
  
  names(gene_pvalues) <- rownames(gene_rank_matrix)
  
  return(gene_pvalues)
  
}