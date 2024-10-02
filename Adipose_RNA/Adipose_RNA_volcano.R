pacman::p_load(tidyverse, readr)

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 18)
)

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

#################
# 12 months
#################

library(readxl)
adipose_rna_limma_output <- read_excel(here::here("../RNAseq/Adipose/crmo12vsalmo12VScrmo0vsalmo0.xlsx")) %>% 
  mutate(`Change in Expression` = if_else(log2FC >0, 'Upregulated', 'Downregulated'),
         FDR_significant = if_else(fdr < 0.05, 'FDR < 0.05', "FDR > 0.05"), 
         Change_and_FDR = if_else(log2FC >0 & fdr < 0.05, "Upregulated", 
                              if_else(log2FC <0 & fdr < 0.05, "Downregulated", "Not Significant")),
         '-log10Pval' = -log10(P.Value), 
         scaled_P  = scale_this(`-log10Pval`), 
         scaled_FC = scale_this(abs(log2FC)), 
         summed_zscores = scaled_P + scaled_FC,
         my_label = if_else(FDR_significant == 'FDR < 0.05', SYMBOL, NA))


adipose_rna_limma_output %>% 
  ggplot(., aes(x = log2FC, y = `-log10Pval`, label = my_label))+
  geom_point(aes(color = Change_and_FDR), alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue2", "gray", 'firebrick1'))+
  ggrepel::geom_label_repel(data = adipose_rna_limma_output %>% 
                              filter(FDR_significant == 'FDR < 0.05') %>% 
                            group_by(Change_and_FDR) %>% 
                              top_n(n = 10, wt = dplyr::desc(summed_zscores)),
                            max.overlaps = Inf,
                              fill = "white")+
  theme(legend.position = 'none')+
  My_Theme+
  xlim(c(-4.3,4.3))


ggsave(here::here("../RNAseq/Adipose/Output/Figures", "CALERIE_12mo_CRvsAL_volcano_CPR.png"))



adipose_rna_limma_output %>% 
  filter(!is.na(my_label))  %>% 
  group_by(Change_and_FDR) %>% 
  dplyr::count()

# # Groups:   Change_and_FDR [2]
# Change_and_FDR     n
# <chr>          <int>
#   1 Downregulated    296
#   2 Upregulated      309

#################
# 24 months
################# 
  adipose_rna_limma_output_24 <- read_excel(here::here("../RNAseq/Adipose/crmo24vsalmo24VScrmo0vsalmo0.xlsx")) %>% 
    mutate(`Change in Expression` = if_else(log2FC >0, 'Upregulated', 'Downregulated'),
           FDR_significant = if_else(fdr < 0.05, 'FDR < 0.05', "FDR > 0.05"), 
           Change_and_FDR = if_else(log2FC >0 & fdr < 0.05, "Upregulated", 
                                    if_else(log2FC <0 & fdr < 0.05, "Downregulated", "Not Significant")),
           '-log10Pval' = -log10(P.Value), 
           scaled_P  = scale_this(P.Value), 
           scaled_FC = scale_this(abs(log2FC)), 
           summed_zscores = scaled_P + scaled_FC,
           my_label = if_else(FDR_significant == 'FDR < 0.05', SYMBOL, NA))
  
  adipose_rna_limma_output_24 %>% 
    ggplot(., aes(x = log2FC, y = `-log10Pval`, label = my_label))+
    geom_point(aes(color = Change_and_FDR), alpha = 0.4)+
    theme_bw()+
    scale_color_manual(values = c("dodgerblue2", "gray", 'firebrick1'))+
    ggrepel::geom_label_repel(fill = "white")+
    theme(legend.position = 'none')+
    My_Theme
  
  ggsave(here::here("../RNAseq/Adipose/Output/Figures", "CALERIE_24mo_CRvsAL_volcano_CPR.png"))
# Plotted 734 - 711 = 23
  
  
  adipose_rna_limma_output_24 %>% filter(!is.na(my_label)) %>% 
  group_by(Change_and_FDR) %>% 
    dplyr::count()
  # Groups:   Change_and_FDR [2]
  # Change_and_FDR     n
  #   <chr>          <int>
  # 1 Downregulated    404
  # 2 Upregulated      330