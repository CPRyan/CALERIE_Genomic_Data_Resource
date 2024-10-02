pacman::p_load(tidyverse, readr, parallel)

numCores <-detectCores()
cl <- makeCluster(numCores)

####################
# Load Data
###################

counts <-read.table(here::here("../RNAseq/Adipose", "age_diet_mrna_counts.3.txt"))
dim(counts)
# counts <-rownames_to_column(counts, "Entrez_Gene_No")
# dim(counts)


pheno_data <-read_table(here::here("../RNAseq/Adipose", "targets.3.txt"))

pheno_data <-pheno_data %>% 
  mutate(Visit.ID = if_else(grepl('BASELINE', FileName),
                            "Baseline", 
                            if_else(grepl('MONTH.12', FileName), 
                                    "12 Month", 
                                    "24 Month")),
         Treatment = if_else(grepl('cr', Subject),
                             "Caloric Restriction",
                             "Ad Libitum"), 
         DEID = sub("_$", "", substr(gsub("X", "", FileName), 1, regexpr("_", FileName) - 1))  
  )





dim(counts)
# 770155    168

# With complete data
counts_nona <-na.omit(counts)
dim(counts_nona)
# 847208    599

nrow(counts)-nrow(counts_nona)
# All na removed anyway.

####################
# All samples
####################

summary_df <-data.frame(rownames(counts), row.names = rownames(counts)); 
names(summary_df)[1] <-"Entrez_Gene_No"

# How many missing per row
summary_df$na_count <- parApply(cl, counts, 1, function(x) sum(is.na(x)))

# mean per row
summary_df$mean <- parApply(cl, counts, 1, mean, na.rm = TRUE)

# sd per row
summary_df$sd <- parApply(cl, counts, 1, sd, na.rm = TRUE)

# median per row
summary_df$median <- parApply(cl, counts, 1, median, na.rm = TRUE)

# min per row
summary_df$min <- parApply(cl, counts, 1, min, na.rm = TRUE)

# max per row
summary_df$max <- parApply(cl, counts, 1, max, na.rm = TRUE)


write_csv(summary_df, here::here("../RNAseq/Adipose/Output", "CALERIE_Adipose_RNAseq_non_na_entrez_gene_list_full_samples.csv"))


#####################
# Baseline Sample
#####################
counts[1:4,1:4]
baseline_barcodes <-pheno_data %>% filter(Visit.ID == 'Baseline') %>% pull(Name)

baseline_counts <-counts[, which(colnames(counts) %in% baseline_barcodes)]

dim(baseline_counts)

baseline_counts_nona <-na.omit(baseline_counts)

dim(baseline_counts_nona)



baseline_summary_df <-data.frame(rownames(baseline_counts), row.names = rownames(baseline_counts)); 
names(baseline_summary_df) <-"Entrez_Gene_No"

# How many missing per row
baseline_summary_df$na_count <- parApply(cl, baseline_counts, 1, function(x) sum(is.na(x)))

# mean per row
baseline_summary_df$mean <- parApply(cl, baseline_counts, 1, mean, na.rm = TRUE)

# sd per row
baseline_summary_df$sd <- parApply(cl, baseline_counts, 1, sd, na.rm = TRUE)

# median per row
baseline_summary_df$median <- parApply(cl, baseline_counts, 1, median, na.rm = TRUE)

# min per row
baseline_summary_df$min <- parApply(cl, baseline_counts, 1, min, na.rm = TRUE)

# max per row
baseline_summary_df$max <- parApply(cl, baseline_counts, 1, max, na.rm = TRUE)

write_csv(baseline_summary_df, here::here("../RNAseq/Adipose/Output", "CALERIE_Adipose_RNAseq_non_na_entrez_gene_list_baseline_samples.csv"))




#####################
# FU1
#####################

counts[1:4,1:4]
fu1_barcodes <-pheno_data %>% filter(Visit.ID == '12 Month') %>% pull(Name)

fu1_counts <-counts[, which(colnames(counts) %in% fu1_barcodes)]


fu1_counts_nona <-na.omit(fu1_counts)

dim(fu1_counts_nona)



fu1_summary_df <-data.frame(rownames(fu1_counts), row.names = rownames(fu1_counts)); 
names(fu1_summary_df) <-"Entrez_Gene_No"

# How many missing per row
fu1_summary_df$na_count <- parApply(cl, fu1_counts, 1, function(x) sum(is.na(x)))

# mean per row
fu1_summary_df$mean <- parApply(cl, fu1_counts, 1, mean, na.rm = TRUE)

# sd per row
fu1_summary_df$sd <- parApply(cl, fu1_counts, 1, sd, na.rm = TRUE)

# median per row
fu1_summary_df$median <- parApply(cl, fu1_counts, 1, median, na.rm = TRUE)

# min per row
fu1_summary_df$min <- parApply(cl, fu1_counts, 1, min, na.rm = TRUE)

# max per row
fu1_summary_df$max <- parApply(cl, fu1_counts, 1, max, na.rm = TRUE)


write_csv(fu1_summary_df, here::here("../RNAseq/Adipose/Output", "CALERIE_Adipose_RNAseq_non_na_entrez_gene_list_12mo_samples.csv"))



#####################
# fu2 Sample
#####################

counts[1:4,1:4]
fu2_barcodes <-pheno_data %>% filter(Visit.ID == '24 Month') %>% pull(Name)

fu2_counts <-counts[, which(colnames(counts) %in% fu2_barcodes)]

fu2_counts_nona <-na.omit(fu2_counts)

dim(fu2_counts_nona)


fu2_summary_df <-data.frame(rownames(fu2_counts), row.names = rownames(fu2_counts)); 
names(fu2_summary_df) <-"Entrez_Gene_No"

# How many missing per row
fu2_summary_df$na_count <- parApply(cl, fu2_counts, 1, function(x) sum(is.na(x)))

# mean per row
fu2_summary_df$mean <- parApply(cl, fu2_counts, 1, mean, na.rm = TRUE)

# sd per row
fu2_summary_df$sd <- parApply(cl, fu2_counts, 1, sd, na.rm = TRUE)

# median per row
fu2_summary_df$median <- parApply(cl, fu2_counts, 1, median, na.rm = TRUE)

# min per row
fu2_summary_df$min <- parApply(cl, fu2_counts, 1, min, na.rm = TRUE)

# max per row
fu2_summary_df$max <- parApply(cl, fu2_counts, 1, max, na.rm = TRUE)


write_csv(fu2_summary_df, here::here("../RNAseq/Adipose/Output", "CALERIE_Adipose_RNAseq_non_na_entrez_gene_list_24mo_samples.csv"))


###############################
