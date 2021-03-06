First we will import the script which merges count data across all
comparisons in all experiments, so we can analyse the raw counts. We’ll
also import the results of the differential expression meta-analysis.

``` r
# source script for merging count data
source('merge_counts.R')
# import differential expression results
ME_results <- read.csv('data/meta_analysis_results.csv', header = T)
```

Now let’s grab some read mapping statistics from the data sets we
imported. For each sample, we’ll extract the number of reads assigned
during the alignment step and the number of reads unassigned.

``` r
library(reshape2, quietly = T)
library(stringr, quietly = T)

# get number of assigned/unassigned reads
assigned_reads <- lapply(all_counts, 
                         function(c) c$stat[c$stat$Status == 'Assigned',-1]
                         ) 
unassigned_reads <- lapply(all_counts, 
                           function(c) apply(c$stat[c$stat$Status != 'Assigned',-1], 
                                             2, sum)
                           ) 

# get metadata columns for our data frame
samples <- lapply(assigned_reads, colnames)
experiments <- lapply(names(all_counts), 
                      function(exp_id) rep(exp_id, ncol(all_counts[[exp_id]]$stat) - 1))

# add all the read mapping and sample information to a single data frame
read_mapping_df <- data.frame('assigned' = unlist(assigned_reads),
                              'unassigned' = unlist(unassigned_reads),
                              'sample' = unlist(samples),
                              'experiment' = unlist(experiments))
read_mapping_df <- melt(read_mapping_df, id.vars = c('sample', 'experiment'),
                        variable.name = 'status', value.name = 'num_reads')

# show a preview of the dataframe
head(read_mapping_df)
```

    ##                                                        sample experiment   status num_reads
    ## 1 X.AV.SS.Data.data.mg989.SRP045441.bam.SRR1544483_1_trim.bam  SRP045441 assigned  27138840
    ## 2 X.AV.SS.Data.data.mg989.SRP045441.bam.SRR1544484_1_trim.bam  SRP045441 assigned  24640229
    ## 3 X.AV.SS.Data.data.mg989.SRP045441.bam.SRR1544485_1_trim.bam  SRP045441 assigned  27593438
    ## 4 X.AV.SS.Data.data.mg989.SRP045441.bam.SRR1544489_1_trim.bam  SRP045441 assigned  29161911
    ## 5 X.AV.SS.Data.data.mg989.SRP045441.bam.SRR1544490_1_trim.bam  SRP045441 assigned  38711564
    ## 6 X.AV.SS.Data.data.mg989.SRP045441.bam.SRR1544491_1_trim.bam  SRP045441 assigned  31663866

Now we can generate a plot showing the percentage of assigned/unassigned
reads for each sample.

``` r
library(ggplot2, quietly = T)

# plot read mapping statistics
read_mapping_plot <- ggplot(read_mapping_df, aes(x = sample, y = num_reads, fill = status)) +
  geom_col(position = 'fill', colour = 'black') +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~experiment, 
             scales = 'free_x',
             space = 'free_x',
             switch = 'x') +
  labs(x = '\nExperiment/Sample', y = 'Percentage of reads\n', 
       title = 'Read mapping statistics, all samples', 
       fill = 'Status') +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 42),
        axis.text.y = element_text(size = 38),
        legend.text = element_text(size = 38),
        legend.title = element_text(size = 42),
        panel.grid = element_blank(),
        plot.title = element_text(size = 50, hjust = 0.5))
        
print(read_mapping_plot)
```

![](plotting_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
# save the figure
ggsave('figures/fig1a.png', read_mapping_plot, width = 40, height = 10)
```

Next we can generate some plots to show similarities between samples. To
visualise sample relationships, we’ll fit a PCA to the full set of
\~16000 genes tested in the meta analysis.

``` r
# import sample metadata
sample_metadata <- read.csv('metadata/RNAseq_sample_metadata.csv', row.names = 1)

# make a gene data frame with different ID types for each gene
human_ids <- org.Hs.eg.db
genedata <- data.frame('ensembl' = count_matrix$gene)
genedata$symbol <- mapIds(human_ids, 
                          genedata$ensembl,
                          column = 'SYMBOL',
                          keytype = 'ENSEMBL')
genedata$entrez <- mapIds(human_ids, 
                          genedata$ensembl,
                          column = 'ENTREZID',
                          keytype = 'ENSEMBL')  

# TMM normalise read counts
fixed_effects <- model.matrix(~ senescence_type, data = sample_metadata)
counts_only <- column_to_rownames(count_matrix, 'gene')
normCountMatrix <- preprocessCounts(counts_only,
                                    sample_metadata,
                                    genedata, 
                                    fixed_effects)

# get log cpms of the counts
cpms <- cpm(normCountMatrix, log = TRUE)

# subset genes that were tested in the meta-analysis
ME_gene_cpms <- cpms[ME_results$ensembl,]

# run PCA on logged counts
countMatrixPC <- prcomp(t(ME_gene_cpms))
countMatrixPC$x <- as.data.frame(countMatrixPC$x)

# add variables for experiment and type of senescence
countMatrixPC$x$experiment <- sample_metadata$experiment
countMatrixPC$x$senescence_type <- sample_metadata$senescence_type

# calculate percentage of variance explained for each PC
countMatrixPC$var_expl <- countMatrixPC$sd^2 / sum(countMatrixPC$sd^2)
```

Now we can plot the results of the PCA, colouring samples by their
experiment ID and type of senescence.

``` r
library(ggplot2, quietly = T)

# scatter plot of samples in PCA space
sample_PCA_plot <- ggplot(countMatrixPC$x, aes(x = PC1, y = PC2, 
                                               colour = experiment, 
                                               shape = senescence_type)) + 
  geom_point(size = 2.5) +   
  labs(x = paste0('PC1 (', round(countMatrixPC$var_expl[1]*100, 2), '%)'), 
       y = paste0('PC2 (', round(countMatrixPC$var_expl[2]*100, 2), '%)'), 
       title = 'First two principal components of TMM-normalised expression values',
       colour = 'Experiment SRP ID',
       shape = 'Type of senescence') +
  theme_bw() +
  theme(text = element_text(size = 18), 
        plot.title = element_text(hjust = 0.5))

print(sample_PCA_plot)
```

![](plotting_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# save the plot
ggsave('figures/pca_plot.png', sample_PCA_plot, width = 12, height = 8)
```

The plot clearly shows some significant batch effects, with samples from
individual experiments clustering together more tightly than samples
with the same type of senescence. To adjust for these effects, we used a
random effects model to model inter-study heterogeneity of log fold
changes.

Next we can extract differentially-expressed genes from the fitted
model. We’ll use a threshold of BH-adjusted p-value \< 0.05 and
\|log2FC\| \> 1 for each type of senescence. Any gene that meets this
threshold for at least one type of senescence will be classified as
“senescence-associated”.

``` r
LFC_threshold <- 1

# subset genes identified as differentially-expressed in oncogene-induced senescence
oncogene_sig_mask <- abs(ME_results$LFC_Oncogene) > LFC_threshold & ME_results$adj_p_Oncogene < 0.05
oncogene_sig_genes <- ME_results[oncogene_sig_mask,]

writeLines(paste(nrow(oncogene_sig_genes), 'genes classified as differentially-expressed in oncogene-induced senescence'))
```

    ## 2768 genes classified as differentially-expressed in oncogene-induced senescence

``` r
# same for replicative senescence 
replicative_sig_mask <- abs(ME_results$LFC_Replicative) > LFC_threshold & ME_results$adj_p_Replicative < 0.05
replicative_sig_genes <- ME_results[replicative_sig_mask,]

writeLines(paste(nrow(replicative_sig_genes), 'genes classified as differentially-expressed in replication-induced senescence'))
```

    ## 1471 genes classified as differentially-expressed in replication-induced senescence

``` r
# same for DNA damage senescence coefficient
DNAdamage_sig_mask <- abs(ME_results$LFC_DNAdamage) > LFC_threshold & ME_results$adj_p_DNAdamage < 0.05
DNAdamage_sig_genes <- ME_results[DNAdamage_sig_mask,]

writeLines(paste(nrow(DNAdamage_sig_genes), 'genes classified as differentially-expressed in DNA damage-induced senescence'))
```

    ## 1588 genes classified as differentially-expressed in DNA damage-induced senescence

``` r
# combine all DE genes into single vector
RNAseq_senescence_genes <- unique(c(oncogene_sig_genes$ensembl, 
                                    replicative_sig_genes$ensembl, 
                                    DNAdamage_sig_genes$ensembl))

writeLines(paste(length(RNAseq_senescence_genes), 'genes classified as differentially-expressed in any form of senescence'))
```

    ## 4091 genes classified as differentially-expressed in any form of senescence

Now we will import the list of genes taken from the innateDB, cellAge,
ageing atlast, and CSgene databases, and join them with the RNA-seq data
to generate a full list of senescence-associated genes.

``` r
library(dplyr, quietly = T)

# import senescence database genes and get LFCs
senescence_db_genes <- read.csv('data/genesets/senescence_db_genes.csv')
senescence_db_genes$entrez <- as.character(senescence_db_genes$entrez)
senescence_db_genes <- left_join(senescence_db_genes, 
                                 genedata[,c('entrez', 'ensembl')],
                                 by = 'entrez')

# filter genes that could not be mapped to ensembl and lookup LFCs
senescence_db_genes <- senescence_db_genes[!is.na(senescence_db_genes$ensembl),]

# combine the full set of senescence genes into a single list for plotting
all_senescence_genes <- list('Oncogene' = oncogene_sig_genes$ensembl,
                             'Replicative' = replicative_sig_genes$ensembl,
                             'DNA-damage' = DNAdamage_sig_genes$ensembl,
                             'Databases' = senescence_db_genes$ensembl)

# get LFCs of senescence database genes
senescence_db_gene_LFCs <- LFClookup(senescence_db_genes$ensembl,
                                     all_senescence_genes[names(all_senescence_genes) != 'Databases'],
                                     column_to_rownames(ME_results, 'ensembl'))
senescence_db_gene_LFCs <- rownames_to_column(senescence_db_gene_LFCs, 'ensembl')
senescence_db_genes <- left_join(senescence_db_genes,
                                 senescence_db_gene_LFCs,
                                 by = 'ensembl')
senescence_db_genes <- column_to_rownames(senescence_db_genes, 'ensembl')

# get LFCs of genes identified via RNA-seq
RNAseq_gene_LFCs <- LFClookup(RNAseq_senescence_genes, 
                              all_senescence_genes, 
                              column_to_rownames(ME_results, 'ensembl'))

RNAseq_gene_LFCs$ensembl <- RNAseq_senescence_genes
senescence_rnaseq_genes <- left_join(RNAseq_gene_LFCs, genedata, by = 'ensembl')
senescence_rnaseq_genes <- column_to_rownames(senescence_rnaseq_genes, 'ensembl')

# merge senescence db genes with rnaseq genes
all_unique_genes <- unique(c(rownames(senescence_rnaseq_genes), 
                             rownames(senescence_db_genes)))
all_genes <- vector(mode = 'list', length = length(all_unique_genes))
for(i in 1:length(all_unique_genes)){
  gene <- all_unique_genes[i]
  if(gene %in% rownames(senescence_db_genes)){
    gene_row <- senescence_db_genes[gene,]
    if(gene %in% rownames(senescence_rnaseq_genes)){
      gene_row$source <- paste(gene_row$source, 'RNAseq', sep = '/')
    }
    all_genes[[i]] <- gene_row
  }
  else if(gene %in% rownames(senescence_rnaseq_genes)){
    gene_row <- senescence_rnaseq_genes[gene,]
    gene_row$source <- 'RNAseq'
    gene_row$senescence_effect <- 'Unknown'
    gene_row$organism <- 'Human'
    all_genes[[i]] <- gene_row
  }
}

# make data frame with all genes
all_genes_df <- do.call(rbind, all_genes)
all_genes_df <- rownames_to_column(all_genes_df, 'ensembl')
all_genes_df$ensembl <- str_replace(all_genes_df$ensembl, '\\..+$', '')

# save the full gene data frame to a csv file
write.csv(all_genes_df, 'data/genesets/all_senescence_genes.csv', row.names = FALSE)
```

To visualise overlap between the different gene sets, we’ll make an
upset plot.

``` r
library(UpSetR, quietly = T)

upset_plot_genesets <- upset(fromList(all_senescence_genes), order.by = 'freq', text.scale = 1.75)

print(upset_plot_genesets)
```

![](plotting_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
# save the figure
png('figures/fig1c.png', width = 2400, height = 2000, units = 'px', res = 300)
upset_plot_genesets
dev.off()
```

    ## png 
    ##   2

Interestingly, the upset plot shows that many of the genes identified as
associated with replicative senescence were also identified as
differentially-expressed in DNA damage-induced senescence. It also shows
that most genes identified as differentially-expressed in senescent
cells are not present in the senescence databases.

We can also make an upset plot showing the overlap between the different
databases.

``` r
senescence_databases <- unique(unlist(str_split(senescence_db_genes$source, '/')))

db_genes_by_db <- lapply(
  senescence_databases, 
  function(db) rownames(senescence_db_genes[grepl(db, senescence_db_genes$source),])
)
names(db_genes_by_db) <- senescence_databases

library(UpSetR, quietly = T)

upset_plot_senescence_dbs <- upset(fromList(db_genes_by_db), order.by = 'freq', text.scale = 1.75)

print(upset_plot_senescence_dbs)
```

![](plotting_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
# save the figure
png('figures/fig1d.png', width = 2400, height = 2000, units = 'px', res = 300)
upset_plot_senescence_dbs
dev.off()
```

    ## png 
    ##   2

Now we can generate some volcano plots using the LFCs estimated from our
mixed model, and their associated (adjusted) p-values.

``` r
library(ggplot2, quietly = TRUE)
library(ggpubr, quietly = TRUE)

# make volcano plots from the results of DE testing
oncogene_vlcplot <- makeVolcanoPlot(ME_results,
                                    'adj_p_Oncogene',
                                    'LFC_Oncogene',
                                    title = 'Oncogene-induced senescence') + 
  xlim(-12, 12) + ylim(0, 20)
  
replicative_vlcplot <- makeVolcanoPlot(ME_results, 
                                       'adj_p_Replicative',
                                       'LFC_Replicative',
                                       title = 'Replication-induced senescence') + 
  xlim(-12, 12) + ylim(0, 20)
DNAdamage_vlcplot <- makeVolcanoPlot(ME_results, 
                                     'adj_p_DNAdamage',
                                     'LFC_DNAdamage', 
                                     title = 'DNA damage-induced senescence') + 
  xlim(-12, 12) + ylim(0, 20)
  
# combine the three volcano plots into a single plot
combined_vlcplot <- ggarrange(oncogene_vlcplot,
                              replicative_vlcplot, 
                              DNAdamage_vlcplot, 
                              nrow = 1, common.legend = T)
                              
print(combined_vlcplot)
```

![](plotting_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# save the figure
ggsave('figures/fig1b.png', combined_vlcplot, width = 18, height = 6)
```
