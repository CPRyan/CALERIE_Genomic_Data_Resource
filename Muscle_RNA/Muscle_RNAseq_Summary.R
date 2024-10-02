pacman::p_load(tidyverse, readr, parallel)

numCores <-detectCores()
cl <- makeCluster(numCores)

####################
# Load Data
###################

counts <-read_csv(here::here("../RNAseq/Muscle", "counts_hg38.csv"))
dim(counts)
names(counts)[1] <-"ensembl_Gene_ID"
# dim(counts)

# No pheno data provided. 
# Try to parse from column names

# First transpose
head(counts)
t_counts <-t(counts[,-1])
t_counts[1:4,1:4]
# Fix names
colnames(t_counts) <-counts$ensembl_Gene_ID
t_counts <-as_tibble(t_counts, rownames = "Bam_ID")
t_counts[1:4,1:4]

# Now parse names
t_counts <-t_counts %>% 
  mutate(Visit.ID = if_else(grepl('Baseline', Bam_ID),
                            "Baseline", 
                            if_else(grepl('12_month', Bam_ID), 
                                    "12 Month", 
                                    "24 Month")),
         DEID = str_extract(Bam_ID,  "(?<=_)(\\d+)(?=_\\w+\\.\\w+)")
  ) %>% 
  select(Bam_ID, Visit.ID, DEID, everything())

head(t_counts)

t_counts %>% group_by(Visit.ID) %>% summarize(n = n())
# # A tibble: 3 × 2
#.   Visit.ID     n
#    <chr>    <int>
#   1 12 Month    44
#   2 24 Month    30
#   3 Baseline    88

# This looks right. But I still have no info on treatment. Annoying I have to go elsewhere for this. 

ids_treat <-haven::read_dta(here::here("../Blood/Data", "CALERIE_DWB220314.dta")) %>% 
  rename(DEID = deidnum) %>% 
  mutate(Treatment = if_else(CR == 1,
                      "Caloric Restriction",
                      "Ad Libitum")) %>% 
  select(DEID, Treatment) %>% 
  sjlabelled::as_factor(DEID) %>% distinct(DEID, .keep_all = TRUE)


new_counts <-inner_join(t_counts, ids_treat, by = "DEID", relationship = 'many-to-many') %>% select(Bam_ID, DEID, Treatment, Visit.ID, everything())

dim(new_counts)
# 162 60609

new_counts[,1:4] %>% group_by(Visit.ID) %>% summarize(n = n())
new_counts[,1:4] %>% group_by(Treatment) %>% summarize(n = n())
new_counts[,1:4] %>% group_by(Treatment, Visit.ID) %>% summarize(n = n())

write_csv(new_counts[,1:4],  here::here("Output", "CALERIE_Muscle_RNAseq_pheno_data_CPR_construction.csv"))


# With complete data
counts_nona <-na.omit(new_counts)
dim(counts_nona)
# 60605   163

nrow(new_counts)-nrow(counts_nona)
# All na removed anyway.

# Make Pheno From transposed counts data
pheno_data <-new_counts[1:4]

# Make nice clean counts data
# Save names
rownamz <-counts$ensembl_Gene_ID
# Convert to matrix.
counts <-as.matrix(counts[,-1])
# Call rownames
rownames(counts)<-rownamz

####################
# All samples
####################

summary_df <-data.frame(rownames(counts), row.names = rownames(counts)); 
names(summary_df)[1] <-"ensemble_Gene_No"

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


write_csv(summary_df, here::here("../RNAseq/Muscle/Output", "CALERIE_Muscle_RNAseq_non_na_ensemble_gene_list_full_samples.csv"))


#####################
# Baseline Sample
#####################
counts[1:4,1:4]
baseline_barcodes <-pheno_data %>% filter(Visit.ID == 'Baseline') %>% pull(Bam_ID)

baseline_counts <-counts[, which(colnames(counts) %in% baseline_barcodes)]

dim(baseline_counts)

baseline_counts_nona <-na.omit(baseline_counts)

dim(baseline_counts_nona)



baseline_summary_df <-data.frame(rownames(baseline_counts), row.names = rownames(baseline_counts)); 
names(baseline_summary_df) <-"ensemble_Gene_No"

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

write_csv(baseline_summary_df, here::here("../RNAseq/Muscle/Output", "CALERIE_Muscle_RNAseq_non_na_ensemble_gene_list_baseline_samples.csv"))




#####################
# FU1
#####################

counts[1:4,1:4]
fu1_barcodes <-pheno_data %>% filter(Visit.ID == '12 Month') %>% pull(Bam_ID)

fu1_counts <-counts[, which(colnames(counts) %in% fu1_barcodes)]


fu1_counts_nona <-na.omit(fu1_counts)

dim(fu1_counts_nona)



fu1_summary_df <-data.frame(rownames(fu1_counts), row.names = rownames(fu1_counts)); 
names(fu1_summary_df) <-"ensemble_Gene_No"

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


write_csv(fu1_summary_df, here::here("../RNAseq/Muscle/Output", "CALERIE_Muscle_RNAseq_non_na_ensemble_gene_list_12mo_samples.csv"))



#####################
# fu2 Sample
#####################

counts[1:4,1:4]
fu2_barcodes <-pheno_data %>% filter(Visit.ID == '24 Month') %>% pull(Bam_ID)

fu2_counts <-counts[, which(colnames(counts) %in% fu2_barcodes)]

fu2_counts_nona <-na.omit(fu2_counts)

dim(fu2_counts_nona)


fu2_summary_df <-data.frame(rownames(fu2_counts), row.names = rownames(fu2_counts)); 
names(fu2_summary_df) <-"ensemble_Gene_No"

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


write_csv(fu2_summary_df, here::here("../RNAseq/Muscle/Output", "CALERIE_Muscle_RNAseq_non_na_ensemble_gene_list_24mo_samples.csv"))


###############################


