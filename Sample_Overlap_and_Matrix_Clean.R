pacman::p_load(tidyverse, readr, UpSetR, sjlabelled, dplyr)


"%nin%" <- function(x, y) {
  return( !(x %in% y) )
}


filter_blfu <- function(df) {
  df %>% 
    group_by(DEID) %>%
    filter(any(Visit.ID == "Baseline") &
             (any(Visit.ID == "12 Month") | any(Visit.ID == "24 Month"))) %>% 
    ungroup()
}


##########################
# Load Data: Adipose DNAm
##########################
pheno_adipose <- read_csv("../Adipose/Output/Data/Corcoran_Pipeline/CALERIE_Adipose_Corcoran_sample_file.csv") %>% 
  mutate(Treatment = if_else(CR == 0, "Ad Libitum", "Caloric Restriction")) %>% 
  select(DEID, Visit.ID, Treatment) %>% na.omit()

##########################
# Load Data: Muscle DNAm
##########################

pheno_muscle <- read_csv("../Muscle/Output/Data/Corcoran_Pipeline/CALERIE_Muscle_Corcoran_sample_file.csv") %>%  
  mutate(Treatment = if_else(CR == 0, "Ad Libitum", "Caloric Restriction")) %>% 
  select(DEID, Visit.ID, Treatment) %>% na.omit()

##########################
# Load Data: Blood DNAm
##########################
read_tsv(here::here("../Blood/Data/phenoData.tsv")) %>% arrange(Sample_Name) %>% pull(Sample_Name)

pheno_blood <- haven::read_dta("../Blood/Data/CALERIE_DWB220314.dta") %>% 
  rename(DEID = deidnum) %>% 
  filter(!is.na(dnamagehannum)) %>% 
  mutate(Visit.ID = if_else(fu == 0, "Baseline", 
                            if_else(fu == 1, '12 Month','24 Month')), 
         Treatment = if_else(CR == 0, "Ad Libitum", "Caloric Restriction")) %>% 
  select(DEID, Visit.ID, Treatment) %>% na.omit()


##########################
# Load Data: Adipose RNAseq
##########################

rna <-read_table(here::here("../RNAseq/Adipose", "targets.3.txt"))

rna <-rna %>% 
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


##########################
# Load Data: Muscle RNAseq
##########################

muscle_rna <-read_csv(here::here("../RNAseq/Muscle/Output", "CALERIE_Muscle_RNAseq_pheno_data_CPR_construction.csv"))


##########################
# Load Data: Muscle small/miRNA
##########################

blood_small <-readxl::read_excel(here::here("../miRNAs", "Plasma_smRNA_sample_description.xlsx"))

##########################
# Load Data: Muscle small/miRNA
##########################

muscle_small <-readxl::read_excel(here::here("../miRNAs", "Muscle_smRNA_sample_description.xlsx"))


##########################
# Load Data: Adipose small/miRNA
##########################

adipose_small <-readxl::read_excel(here::here("../miRNAs", "Adipose_smRNA_sample_description.xlsx"))


##########################
# Load Data: SNPs
##########################

snps <-haven::read_dta(here::here("../SNPs", "CALERIE_GeneticPCs200727.dta")) %>% 
  rename(DEID = deidnum) 


##############################################################################
# Filter to data with baseline and at least one follow-up
# This step is no longer used in the final MS, but is provided in the Supplementary Materials.
##############################################################################

# snps 
# not necessary

# pheno_blood2 <-filter_blfu(pheno_blood) %>% filter(DEID != "10015" & DEID != "24270")
# pheno_adipose2 <-filter_blfu(pheno_adipose) %>% filter(DEID != "10015" & DEID != "24270")
# pheno_muscle2 <-filter_blfu(pheno_muscle) %>% filter(DEID != "10015" & DEID != "24270")
# rna2 <-filter_blfu(rna) %>% filter(DEID != "10015" & DEID != "24270")
# muscle_rna2 <-filter_blfu(muscle_rna) %>% filter(DEID != "10015" & DEID != "24270")
# muscle_small2 <-filter_blfu(muscle_small) %>% filter(DEID != "10015" & DEID != "24270")
# adipose_small2 <-filter_blfu(adipose_small) %>% filter(DEID != "10015" & DEID != "24270")
# blood_small2 <-filter_blfu(blood_small) %>% filter(DEID != "10015" & DEID != "24270")


#####################################
# Filter off 10015 and 24270 - withdrew consent
#####################################

pheno_blood2 <-pheno_blood %>% filter(DEID != "10015" & DEID != "24270")
pheno_adipose2 <-pheno_adipose %>% filter(DEID != "10015" & DEID != "24270")
pheno_muscle2 <-pheno_muscle %>% filter(DEID != "10015" & DEID != "24270")
rna2 <-rna %>% filter(DEID != "10015" & DEID != "24270")
muscle_rna2 <-muscle_rna %>% filter(DEID != "10015" & DEID != "24270")
muscle_small2 <-muscle_small %>% filter(DEID != "10015" & DEID != "24270")
adipose_small2 <-adipose_small %>% filter(DEID != "10015" & DEID != "24270")
blood_small2 <-blood_small %>% filter(DEID != "10015" & DEID != "24270")

##########################
# Subset by treatment for upsetplot
##########################

adiposeCR <-pheno_adipose2 %>% filter(Treatment == "Caloric Restriction")
adiposeAL <-pheno_adipose2 %>% filter(Treatment == "Ad Libitum")

muscleCR <-pheno_muscle2 %>% filter(Treatment == "Caloric Restriction")
muscleAL <-pheno_muscle2 %>% filter(Treatment == "Ad Libitum")

bloodCR <-pheno_blood2 %>% filter(Treatment == "Caloric Restriction")
bloodAL <-pheno_blood2 %>% filter(Treatment == "Ad Libitum")

rnaCR <-rna2 %>% filter(Treatment == "Caloric Restriction")
rnaAL <-rna2 %>% filter(Treatment == "Ad Libitum")

muscle_rnaCR <-muscle_rna2 %>% filter(Treatment == "Caloric Restriction")
muscle_rnaAL <-muscle_rna2 %>% filter(Treatment == "Ad Libitum")

blood_smallCR <-blood_small2 %>% filter(TREATMENT == "Caloric Restriction")
blood_smallAL <-blood_small2 %>% filter(TREATMENT == "Ad Libitum")

muscle_smallCR <-muscle_small2 %>% filter(TREATMENT == "Caloric Restriction")
muscle_smallAL <-muscle_small2 %>% filter(TREATMENT == "Ad Libitum")

adipose_smallCR <-adipose_small2 %>% filter(TREATMENT == "Caloric Restriction")
adipose_smallAL <-adipose_small2 %>% filter(TREATMENT == "Ad Libitum")


SNPs <-snps %>% 
  filter(DEID %in% c(adiposeCR$DEID, muscleCR$DEID, bloodCR$DEID, rnaCR$DEID, muscle_rnaCR$DEID, blood_smallCR$DEID, muscle_smallCR$DEID, adipose_smallCR$DEID, 
                      adiposeAL$DEID, muscleAL$DEID, bloodAL$DEID, rnaAL$DEID, muscle_rnaAL$DEID, blood_smallAL$DEID, muscle_smallAL$DEID, adipose_smallAL$DEID)) %>% pull(DEID)


my_list <-list('Adipose CR DNAm' = adiposeCR$DEID, 
               'Adipose AL DNAm' = adiposeAL$DEID,
               "Muscle CR DNAm" = muscleCR$DEID,  
               "Muscle AL DNAm" = muscleAL$DEID, 
               "Blood CR DNAm" = bloodCR$DEID, 
               "Blood AL DNAm" = bloodAL$DEID, 
               "Adipose CR mRNA" = rnaCR$DEID, 
               "Adipose AL mRNA" = rnaAL$DEID, 
               "Muscle CR mRNA" = muscle_rnaCR$DEID,
               "Muscle AL mRNA" = muscle_rnaAL$DEID,
               "Blood CR smRNAs" = blood_smallCR$DEID, 
               "Blood AL smRNAs" = blood_smallAL$DEID, 
               "Muscle CR smRNAs" = muscle_smallCR$DEID, 
               "Muscle AL smRNAs" = muscle_smallAL$DEID, 
               "Adipose CR smRNAs" = adipose_smallCR$DEID, 
               "Adipose AL smRNAs" = adipose_smallAL$DEID, 
               "SNPs" = SNPs[SNPs != "86405"]) # *** See below for why removed (essentially missing Visit.ID and Treatment.)




upsetplot <-upset(fromList(my_list), 
                  keep.order = TRUE,
                  order.by = c('freq'),
                  main.bar.color = 'navy', 
                  text.scale = 2.2, 
                  set_size.show = TRUE,
                  mb.ratio = c(0.4, 0.6),
                  set_size.scale_max = 250,
                  sets = c('Adipose AL mRNA', 'Adipose AL smRNAs', 'Adipose AL DNAm', 
                           'Muscle AL mRNA', 'Muscle AL smRNAs', 'Muscle AL DNAm', 
                           'Blood AL smRNAs', 'Blood AL DNAm', 
                           'Adipose CR mRNA','Adipose CR smRNAs', 'Adipose CR DNAm', 
                           'Muscle CR mRNA', 'Muscle CR smRNAs', 'Muscle CR DNAm', 
                           'Blood CR smRNAs', 'Blood CR DNAm', 
                           'SNPs'),
                  sets.bar.color=c('orange1','orange3', 'orange4',
                                   "cadetblue1","cadetblue3",'cadetblue4',
                                   "red1","red4",
                                   'orange1','orange3', 'orange4',
                                   "cadetblue1","cadetblue3",'cadetblue4',
                                   "red1","red4",
                                   'magenta4'))


png(here::here("Output/Figures", "CALERIE_DNAm_mRNA_smRNA_SNPs_BLplusOneFU_minus10015_24270_or_snps_only.png"), units = "in", height = 7, width = 13, res = 300)
print(upsetplot)                               ## use print here
dev.off()


metadata <-data.frame(name = c('Adipose AL mRNA', 'Adipose AL smRNAs', 'Adipose AL DNAm', 
                       'Muscle AL mRNA', 'Muscle AL smRNAs', 'Muscle AL DNAm', 
                       'Blood AL smRNAs', 'Blood AL DNAm', 
                       'Adipose CR mRNA','Adipose CR smRNAs', 'Adipose CR DNAm', 
                       'Muscle CR mRNA', 'Muscle CR smRNAs', 'Muscle CR DNAm', 
                       'Blood CR smRNAs', 'Blood CR DNAm', 
                       'SNPs'), 
                     group =  c("AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", 
                        "CR", "CR", "CR", "CR", "CR", "CR", "CR", "CR", "Both"))
  
upsetplot <-upset(fromList(my_list), 
                  keep.order = TRUE,
                  order.by = c('freq'),
                  main.bar.color = 'navy', 
                  text.scale = 2.2, 
                  set_size.show = TRUE,
                  mb.ratio = c(0.4, 0.6),
                  set_size.scale_max = 250,
                  sets = c('Adipose AL mRNA', 'Adipose CR mRNA',
                           'Adipose AL smRNAs','Adipose CR smRNAs', 
                           'Adipose AL DNAm', 'Adipose CR DNAm', 
                           'Muscle AL mRNA', 'Muscle CR mRNA', 
                           'Muscle AL smRNAs', 'Muscle CR smRNAs', 
                           'Muscle AL DNAm', 'Muscle CR DNAm', 
                           'Blood AL smRNAs', 'Blood CR smRNAs', 
                           'Blood AL DNAm', 'Blood CR DNAm', 
                           'SNPs'),
                  sets.bar.color=c('orange1','orange1',
                                   'orange3', 'orange3',
                                   'orange4', 'orange4',
                                   "cadetblue1","cadetblue1",
                                   "cadetblue3","cadetblue3",
                                   'cadetblue4','cadetblue4',
                                   "red4","red4",
                                   "red1", "red1",
                                   'magenta4'), 
                  set.metadata = list(data = metadata, plots = list(
                    list(
                      type = "matrix_rows",
                      column = "group",
                      colors = c(AL = "navy", CR = "red", Both = "gray")))))

png(here::here("Output/Figures", "CALERIE_DNAm_mRNA_smRNA_SNPs_BLplusOneFU_minus10015_24270_or_snps_only_tissue_ordered_groupedcolor.png"), units = "in", height = 7, width = 13, res = 300)
print(upsetplot)                               ## use print here
dev.off()



##########################
# Create Full Cross Matrix
##########################

# Pull of all names from all datasets (to make sure I don't miss any)
all_ids <-bind_rows(pheno_blood %>% filter(DEID != "24270") %>% as_character(DEID), 
                    pheno_adipose  %>% as_character(DEID), 
                    pheno_muscle  %>% as_character(DEID), 
                    rna  %>% as_character(DEID), 
                    muscle_rna %>% as_character(DEID),
                    muscle_small %>% as_character(DEID),
                    adipose_small %>% as_character(DEID),
                    blood_small %>% filter(DEID != "10015" & DEID != "24270") %>% as_character(DEID),
                    snps %>% as_character(DEID) %>% filter(DEID !="86405")) %>% # Has no other molecular data
  select(DEID, Visit.ID, Treatment) %>% 
  distinct(DEID, .keep_all = TRUE)


SNPs2 <-snps %>% filter(DEID !="86405") %>% pull(DEID)

# 218 (this includes the 2 samples with data only for blood small RNAs, and two samples that withdrew consent (10015 and 24270))
all_ids <-all_ids %>% filter(DEID != "10015" & DEID != "24270")
# 218

all_list <-list("all_IDs" = all_ids$DEID,
                "Blood baseline DNAm" = pheno_blood %>% filter(Visit.ID == "Baseline") %>% filter(DEID != "24270") %>% pull(DEID), 
                "Blood 12mo DNAm" = pheno_blood %>% filter(Visit.ID == "12 Month") %>% pull(DEID), 
                "Blood 24mo DNAm" = pheno_blood %>% filter(Visit.ID == "24 Month") %>% pull(DEID), 
                'Adipose baseline DNAm' = pheno_adipose %>% filter(Visit.ID == "Baseline") %>% pull(DEID), 
                'Adipose 12mo DNAm' = pheno_adipose %>% filter(Visit.ID == "12 Month") %>% pull(DEID),
                'Adipose 24mo DNAm' = pheno_adipose %>% filter(Visit.ID == "24 Month") %>% pull(DEID),
                "Muscle baseline DNAm" = pheno_muscle %>% filter(Visit.ID == "Baseline") %>% pull(DEID),  
                "Muscle 12mo DNAm" = pheno_muscle %>% filter(Visit.ID == "12 Month") %>% pull(DEID), 
                "Muscle 24mo DNAm" = pheno_muscle %>% filter(Visit.ID == "24 Month") %>% pull(DEID), 
                "Adipose baseline mRNA" = rna %>% filter(Visit.ID == "Baseline") %>% pull(DEID), 
                "Adipose 12mo mRNA" = rna %>% filter(Visit.ID == "12 Month") %>% pull(DEID), 
                "Adipose 24mo mRNA" = rna %>% filter(Visit.ID == "24 Month") %>% pull(DEID),
                "Muscle baseline mRNA" = muscle_rna %>% filter(Visit.ID == "Baseline") %>% pull(DEID), 
                "Muscle 12mo mRNA" = muscle_rna %>% filter(Visit.ID == "12 Month") %>% pull(DEID), 
                "Muscle 24mo mRNA" = muscle_rna %>% filter(Visit.ID == "24 Month") %>% pull(DEID),
                "Blood baseline small RNAs" = blood_small %>% filter(Visit.ID == "Baseline") %>% filter(DEID != "10015" & DEID != "24270") %>% pull(DEID), 
                "Blood 12mo small RNAs" = blood_small %>% filter(Visit.ID == "12 Month") %>% pull(DEID), 
                "Blood 24mo small RNAs" = blood_small %>% filter(Visit.ID == "24 Month") %>% pull(DEID), 
                'Adipose baseline small RNAs' = adipose_small %>% filter(Visit.ID == "Baseline") %>% pull(DEID), 
                'Adipose 12mo small RNAs' = adipose_small %>% filter(Visit.ID == "12 Month") %>% pull(DEID),
                'Adipose 24mo small RNAs' = adipose_small %>% filter(Visit.ID == "24 Month") %>% pull(DEID),
                "Muscle baseline small RNAs" = muscle_small %>% filter(Visit.ID == "Baseline") %>% pull(DEID),  
                "Muscle 12mo small RNAs" = muscle_small %>% filter(Visit.ID == "12 Month") %>% pull(DEID), 
                "Muscle 24mo small RNAs" = muscle_small %>% filter(Visit.ID == "24 Month") %>% pull(DEID), 
                "SNPs" = SNPs2[SNPs2 != "86405" & SNPs2 != "24270"])

matrix_list <-fromList(all_list)

dim(matrix_list)

# First bind up what I have into a matrix with phenotypic info
bound_matrix <-bind_cols(all_ids, matrix_list) %>% 
  # Remove anyone missing visit.id or treatment id
  filter(!is.na(Visit.ID) | !is.na(Treatment)) 


# Flatten and search for anyone in the list that aren't in the bound matrix. (i.e. anyone who doesn't have Visit.ID or Treatment)
flatten(all_list)[flatten(all_list) %nin% bound_matrix$DEID]


# *** See above - removed from final.
bind_cols(all_ids, matrix_list)  %>% filter(DEID == "86405")
# Went back above and removed. 86405 not present now.



# Both of the two blood small RNAs do not have treatment group. This needs to be fixed. 
# Find the IDS with no treatment info
bound_matrix %>% filter(is.na(Treatment)) %>% View()

# 48097
# 82411

# Both are blood small RNAs data
# Checked records Mega_merge_uber_pheno_ordered.txt sent by Melissa on April 11th. 
# 48097 == AL
# 82411 == CR

bound_matrix[bound_matrix$DEID == '48097',]$Treatment <- 'Ad Libitum'
bound_matrix[bound_matrix$DEID == '82411',]$Treatment <- 'Caloric Restriction'
# Fixed. Do not worry about the upset plot above, it has treatment correctly.

























# Table of Matrix results  
bound_matrix %>% 
  group_by(Treatment, Visit.ID) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  filter(Treatment == "Caloric Restriction") %>% 
  janitor::adorn_totals("row") %>% 
  filter(Treatment == "Total") %>% View()

# Table of Matrix results  
bound_matrix %>% 
  group_by(Treatment, Visit.ID) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  filter(Treatment == "Ad Libitum") %>% 
  janitor::adorn_totals("row") %>% 
  filter(Treatment == "Total") %>% View()





bound_matrix %>% 
  select(-Visit.ID, -all_IDs) %>% 
  write_csv(here::here("..", "CALERIE_Molecular_Sample_Matrix_2024_08_01.csv"))

# A note about the upset plot and the table of these results:
# These don't necessarily add up in some places where it might seem like they should. 
# For example, the left hand side of the plot shows 143 for Blood CR DNam. The number for baseline blood for CR is 142. 
# That's fine! 
# 143 is the total number. One person had data for 12mo and 24 mo, but not baseline. Similar things exist elsewhere. 
# ALSO
# The one sample I kick out that is missing data for all variables 86405
# THIS IS FINE. 
# Seems as though someone in the trail did not consent to provide molecular data. So ultimately the N in the trail (218) does not == the number with molecular data (217)
# 