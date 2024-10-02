
library(haven)
library(tidyverse)
library(sjlabelled)
library(labelled)
library(gage)
library(fgsea)


My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 18)
)

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


source(here::here("Code/Other", "GSEA_Function.R"))
source(here::here("Code/Other", "plot_genesets_clusters_Function.R"))
# Take Adipose RNA
# Rank by score (delta beta and pval std)
# Feed into GSEA function



GO_file = here::here("Data/Other", "c5.go.bp.v7.5.symbols.gmt")



library(readxl)
adipose_rna_limma_output <- read_excel(here::here("Data/Other/Adipose_RNAseq/crmo12vsalmo12VScrmo0vsalmo0.xlsx")) %>% 
  mutate(`Change in Expression` = if_else(log2FC >0, 'Upregulated', 'Downregulated'),
         FDR_significant = if_else(fdr < 0.05, 'FDR < 0.05', "FDR > 0.05"), 
         Change_and_FDR = if_else(log2FC >0 & fdr < 0.05, "Upregulated", 
                                  if_else(log2FC <0 & fdr < 0.05, "Downregulated", "Not Significant")),
         '-log10Pval' = -log10(P.Value), 
         scaled_P  = scale_this(`-log10Pval`), 
         scaled_FC = scale_this(abs(log2FC)), 
         summed_zscores = scaled_P + scaled_FC,
         my_label = if_else(FDR_significant == 'FDR < 0.05', SYMBOL, 'NA'))

################################
# Make gene list
################################

gene_list <-adipose_rna_limma_output$log2FC
names(gene_list) <-adipose_rna_limma_output$SYMBOL

gene_list <-sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

head(gene_list)
tail(gene_list)

################################
# Run GSEA
################################
res = GSEA(gene_list, GO_file, pval = 0.05)
res$Plot

dim(res$Results)


res$Plot

saveRDS(res$Results, here::here("Output/Data/Other/Enrichment/CALERIE_adipose_rnaseq_12month_enrichment_results.rds"))
ggsave(here::here("Output/Figures/Enrichment/Adipose_RNAseq", "CALERIE_adipose_rnaseq_12month_enrichment.png"), width = 6, height = 5)
################################
# Run geneset clusters
################################
plot_geneset_clusters( gs_results = res$Results[res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )

plot_geneset_clusters( gs_results = res$Results[res$Results$NES < 0, ], 
                       main = "Down-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )

################################
# 24 months
################################
adipose_rna_limma_output_24 <- read_excel(here::here("Data/Other/Adipose_RNAseq/crmo24vsalmo24VScrmo0vsalmo0.xlsx")) %>% 
  mutate(`Change in Expression` = if_else(log2FC >0, 'Upregulated', 'Downregulated'),
         FDR_significant = if_else(fdr < 0.05, 'FDR < 0.05', "FDR > 0.05"), 
         Change_and_FDR = if_else(log2FC >0 & fdr < 0.05, "Upregulated", 
                                  if_else(log2FC <0 & fdr < 0.05, "Downregulated", "Not Significant")),
         '-log10Pval' = -log10(P.Value), 
         scaled_P  = scale_this(`-log10Pval`), 
         scaled_FC = scale_this(abs(log2FC)), 
         summed_zscores = scaled_P + scaled_FC,
         my_label = if_else(FDR_significant == 'FDR < 0.05', SYMBOL, 'NA'))

################################
# Make gene list
################################

gene_list <-adipose_rna_limma_output_24$log2FC
names(gene_list) <-adipose_rna_limma_output_24$SYMBOL

gene_list <-sort(gene_list, decreasing = FALSE)
gene_list = gene_list[!duplicated(names(gene_list))]

head(gene_list)
tail(gene_list)

################################
# Run GSEA
################################
res = GSEA(gene_list, GO_file, pval = 0.05)

res <-res %>% mutate(pathway)

dim(res$Results)

res$Plot
saveRDS(res$Results, here::here("Output/Data/Other/Enrichment/CALERIE_adipose_rnaseq_24month_enrichment_results.rds"))
ggsave(here::here("Output/Figures/Enrichment/Adipose_RNAseq", "CALERIE_adipose_rnaseq_24month_enrichment.png"), width = 6, height = 5)

################################
# Run geneset clusters
################################
plot_geneset_clusters( gs_results = res$Results[res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )

plot_geneset_clusters( gs_results = res$Results[res$Results$NES < 0, ], 
                       main = "Down-regulated GSEA clusters",
                       GO_file = GO_file,
                       min.sz = 4 )