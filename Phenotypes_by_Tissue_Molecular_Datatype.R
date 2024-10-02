#######
# Load Packages
#######
pacman::p_load(tidyverse, haven, sjlabelled, rlang)


#######
# Load data
#######
# Load sample availability matrix
mat <-read_csv(here::here("CALERIE_Molecular_Sample_Matrix_2024_04_23.csv")) 

sum(mat[,-c(1:2)])
sum(colSums(mat[,-c(1:2)]))


# Load phenotype Matrix
Mega_merge_uber_pheno_ordered <- read_delim("Mega_merge_uber_pheno_ordered.txt",
                                            delim = "\t", escape_double = FALSE,
                                            ) %>% 
  filter(Metabolon.Timepoint == "Baseline") %>% 
  select(DEIDNUM, AGEBL, GENDER, race3, mbmi)


# This fixes some missing data for 65202
dan <-read_dta(here::here("../../RMR_analysis/Data", "CALERIE_DWB220314.dta"))
fix65202 <-dan %>% filter(deidnum == "65202") %>% select(deidnum, sex, bage, mbmi, contains('race'), ethnic) %>% slice(1)



########
# Merge and pull variables
########

df <-left_join(mat, Mega_merge_uber_pheno_ordered, by = c("DEID" = "DEIDNUM")) 

df[df$DEID == "65202", "AGEBL"] <-fix65202$bage
df[df$DEID == "65202", "GENDER"] <-"Female"
df[df$DEID == "65202", "race3"] <-"White"
df[df$DEID == "65202", "mbmi"] <-fix65202$mbmi

# Function to filter and summarize data based on complex conditions
summarize_data <- function(data, condition) {
  data %>%
    filter(!!parse_expr(condition)) %>%  # Apply the filtering condition
    reframe(
      N = n(), # Number of rows meeting the condition
      mean_AGEBL = mean(AGEBL),
      sd_AGEBL = sd(AGEBL), 
      range_AGEBL = range(AGEBL),
      percent_female = mean(GENDER == "Female") * 100,
      percent_white = mean(race3 == "White") * 100,
      percent_black = mean(race3 == "Black") * 100,
      percent_other = mean(race3 == "Other") * 100,
      mean_mbmi = mean(mbmi),
      sd_mbmi = sd(mbmi)
    )
}



# Any sample (baseline, 12mo, 24mo)
summarize_data(df, "SNPs == 1")
summarize_data(df, "`Blood baseline DNAm` == 1 | `Blood 12mo DNAm` == 1 | `Blood 24mo DNAm` == 1")
summarize_data(df, "`Blood baseline small RNAs` == 1 | `Blood 12mo small RNAs` == 1 | `Blood 24mo small RNAs` == 1")
summarize_data(df, "`Adipose baseline DNAm` == 1 | `Adipose 12mo DNAm` == 1 | `Adipose 24mo DNAm` == 1")
summarize_data(df, "`Adipose baseline mRNA` == 1 | `Adipose 12mo mRNA` == 1 | `Adipose 24mo mRNA` == 1")
summarize_data(df, "`Adipose baseline small RNAs` == 1 | `Adipose 12mo small RNAs` == 1 | `Adipose 24mo small RNAs` == 1")
summarize_data(df, "`Muscle baseline DNAm` == 1 | `Muscle 12mo DNAm` == 1 | `Muscle 24mo DNAm` == 1")
summarize_data(df, "`Muscle baseline mRNA` == 1 | `Muscle 12mo mRNA` == 1 | `Muscle 24mo mRNA` == 1")
summarize_data(df, "`Muscle baseline small RNAs` == 1 | `Muscle 12mo small RNAs` == 1 | `Muscle 24mo small RNAs` == 1")


# Baseline and at least one follow-up (12mo or 24mo)
summarize_data(df, "SNPs == 1")
summarize_data(df, "`Blood baseline DNAm` == 1 & `Blood 12mo DNAm` == 1 | `Blood baseline DNAm` == 1 & `Blood 24mo DNAm` == 1")
summarize_data(df, "`Blood baseline small RNAs` == 1 & `Blood 12mo small RNAs` == 1 | `Blood baseline small RNAs` == 1 & `Blood 24mo small RNAs` == 1")
summarize_data(df, "`Adipose baseline DNAm` == 1 & `Adipose 12mo DNAm` == 1 | `Adipose baseline DNAm` == 1 & `Adipose 24mo DNAm` == 1")
summarize_data(df, "`Adipose baseline mRNA` == 1 & `Adipose 12mo mRNA` == 1 | `Adipose baseline mRNA` == 1 & `Adipose 24mo mRNA` == 1")
summarize_data(df, "`Adipose baseline small RNAs` == 1 & `Adipose 12mo small RNAs` == 1 | `Adipose baseline small RNAs` == 1 & `Adipose 24mo small RNAs` == 1")
summarize_data(df, "`Muscle baseline DNAm` == 1 & `Muscle 12mo DNAm` == 1 | `Muscle baseline DNAm` == 1 & `Muscle 24mo DNAm` == 1")
summarize_data(df, "`Muscle baseline mRNA` == 1 & `Muscle 12mo mRNA` == 1 | `Muscle baseline mRNA` == 1 & `Muscle 24mo mRNA` == 1")
summarize_data(df, "`Muscle baseline small RNAs` == 1 & `Muscle 12mo small RNAs` == 1 | `Muscle baseline small RNAs` == 1 & `Muscle 24mo small RNAs` == 1")

