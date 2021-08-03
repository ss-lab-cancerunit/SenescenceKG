library(ChIPpeakAnno)
library(biomaRt)
library(dplyr)
library(stringr)

promoters <- read.table('datasets/TFs/allPromoters_roadmap.bed', header = FALSE, sep = '\t',
                        col.names = c('Chr', 'Start', 'End'))
promoters <- group_split(promoters, Chr) %>%
  lapply(function(df){
    df$name <- paste(df$Chr, 1:nrow(df), sep = '_')
    return(df)}) %>%
  do.call(rbind, .)

promoters_grange <- toGRanges(promoters, header = TRUE)

hg_mart <- useMart('ensembl', 'hsapiens_gene_ensembl')

promoter_annotations <- annotatePeakInBatch(promoters_grange, hg_mart, featureType = 'TSS')

promoters_with_genes <- as.data.frame(promoter_annotations@elementMetadata)
promoters_with_genes$chr <- as.vector(promoter_annotations@seqnames)

write.table(promoters_with_genes, file = 'datasets/TFs/gene_promoters.tsv', sep = '\t',
            row.names = FALSE)