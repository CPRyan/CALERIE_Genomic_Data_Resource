library(haven)
library(tidyverse)
library(sjPlot)
library(GGally)

pcclocks <- read_csv("../../Nature_Aging_Figures/Data/CALERIE_PC-clocks_agedwb_CPR_fixed.csv") %>% 
  select(barcode, PCHorvath1, PCHorvath1Resid, PCHannum, PCHannumResid, PCPhenoAge, PCPhenoAgeResid, PCGrimAge, PCGrimAgeResid)
clocks <-read_csv("../../Nature_Aging_Figures/Data/CALERIE_agedwb_horvath_datout.output.csv") %>% 
  select(barcode, Female, Age, CR, fu, sample_name, DNAmAge, AgeAccelerationResidual, DNAmAgeHannum, AgeAccelerationResidualHannum, DNAmPhenoAge, AgeAccelPheno, DNAmGrimAge, AgeAccelGrim)
pace <-read_table("../../Nature_Aging_Figures/Data/CALERIE_PoAm45.txt") %>% rename(DunedinPACE = PoAm45) %>% select(ID, DunedinPACE)

calerie <-left_join(left_join(clocks, pcclocks, by = "barcode"), pace, by=c("sample_name" = "ID") ) %>% 
  filter(fu == 'base') %>% 
  mutate(DunedinPACEResid = residuals(lm(DunedinPACE ~ Age))) %>% 
  na.omit()

my_cols <- c("navy", "red")  

# First, residualized clocks
source(here::here("Dayoon_Scatterplot_Function.R"))

# Only PC clocks for main manuscript
agevar = c("Age", "PCHorvath1Resid", "PCHannumResid", "PCPhenoAgeResid", "PCGrimAgeResid", "DunedinPACEResid")

label = c("PCHorvath1Resid" = "PCHorvath1Resid", 
          "PCHannumResid" = "PCHannumResid", 
          "PCPhenoAgeResid" = "PCPhenoAgeResid", 
          "PCGrimAgeResid" = "PCGrimAgeResid", 
          "DunedinPACEResid" = "DunedinPACEResid") 

axis_type = c("PCHorvath1Resid"= "float", "PCHannumResid"= "float", "PCPhenoAgeResid"= "float", "PCGrimAgeResid"= "float", "DunedinPACEResid" = "float")

pdf(here::here("Output/Figures", "Scatterplot_Matrix_Clocks_Age_lm_slope_RESIDUALS_for_Biorepository.pdf"), width = 10, height = 10)
plot_baa(data = calerie, agevar = agevar, label = label, axis_type = axis_type)
dev.off()



# All clocks for supplements
agevar = c("AgeAccelerationResidual", "PCHorvath1Resid", "AgeAccelerationResidualHannum", "PCHannumResid", "AgeAccelPheno", "PCPhenoAgeResid", "AgeAccelGrim", "PCGrimAgeResid", "DunedinPACEResid")

label = c("AgeAccelerationResidual" = "HorvathResid",
          "PCHorvath1Resid" = "PCHorvath1Resid", 
          "AgeAccelerationResidualHannum" = "HannumResid",
          "PCHannumResid" = "PCHannumResid", 
          "AgeAccelPheno" = "PhenoAgeResid",
          "PCPhenoAgeResid" = "PCPhenoAgeResid",  
          "AgeAccelGrim" = "GrimAgeResid",
          "PCGrimAgeResid" = "PCGrimAgeResid", 
          "DunedinPACEResid" = "DunedinPACEResid") 

axis_type = c("AgeAccelerationResidual" = 'float', "PCHorvath1Resid"= "float", "AgeAccelerationResidualHannum" = 'float', "PCHannumResid"= "float", "AgeAccelPheno"  = 'float', "PCPhenoAgeResid"= "float", "AgeAccelGrim"  = 'float', "PCGrimAgeResid"= "float", "DunedinPACEResid" = "float")

pdf(here::here("Output/Figures", "Scatterplot_Matrix_Clocks_Age_lm_slope_RESIDUALS_ALL_CLOCKS_for_Biorepository.pdf"), width = 15, height = 15)
plot_baa(data = calerie, agevar = agevar, label = label, axis_type = axis_type)
dev.off()

#####################
# Now age correlations
#####################
core_clocks <-calerie %>% 
  select(Age, CR, PCHorvath1, PCHannum, PCPhenoAge, PCGrimAge, DunedinPACE)

all_clocks <-calerie %>% 
  select(Age, CR, DNAmAge, PCHorvath1, DNAmAgeHannum, PCHannum, DNAmPhenoAge, PCPhenoAge, DNAmGrimAge, PCGrimAge, DunedinPACE)

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 12))

core_clocks %>% 
  pivot_longer(-c(Age, CR), names_to = "Clock_Type", values_to = "Clock_Value") %>% 
  mutate(Clock_Type = forcats::fct_relevel(Clock_Type, "PCHorvath1", "PCHannum", "PCPhenoAge", "PCGrimAge", "DunedinPACE")) %>% 
  ggplot(., aes(x = Age, y = `Clock_Value`, color = `Clock_Type`))+
  geom_point(aes(fill = Clock_Type), pch = 21, color = 'black', alpha = 0.6)+
  geom_abline(aes(intercept = 0, slope = 1, color = Clock_Type), linetype = 'dashed')+
  theme_bw()+
  facet_wrap(~`Clock_Type`, nrow = 5, scales = 'free')+
  facet_grid(Clock_Type ~ ., scales = 'free', switch = "x") +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "inside",) +
  ylab("Epigenetic Age/Pace of Aging")+
  paletteer::scale_color_paletteer_d("fishualize::Etheostoma_spectabile")+
  paletteer::scale_fill_paletteer_d("fishualize::Etheostoma_spectabile")+
  My_Theme

ggsave(here::here("Output/Figures", "CALERIE_clocks_by_age.png"), width = 3, height = 9)



all_clocks %>% 
  rename(Horvath = DNAmAge, 
         Hannum = DNAmAgeHannum, 
         PhenoAge = DNAmPhenoAge, 
         GrimAge = DNAmGrimAge) %>% 
  pivot_longer(-c(Age, CR), names_to = "Clock_Type", values_to = "Clock_Value") %>% 
  mutate(Clock_Type = forcats::fct_relevel(Clock_Type, "Horvath", "PCHorvath1","Hannum", "PCHannum", "PhenoAge", "PCPhenoAge", "GrimAge", "PCGrimAge", "DunedinPACE")) %>% 
  ggplot(., aes(x = Age, y = `Clock_Value`, color = `Clock_Type`))+
  geom_point(aes(fill = Clock_Type), pch = 21, color = 'black', alpha = 0.6)+
  geom_abline(aes(intercept = 0, slope = 1, color = Clock_Type), linetype = 'dashed')+
  theme_bw()+
  facet_wrap(~`Clock_Type`, scales = 'free')+
  facet_wrap(Clock_Type ~ ., scales = 'free', switch = "x") +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "inside",) +
  ylab("Epigenetic Age/Pace of Aging")+
  paletteer::scale_color_paletteer_d("ggthemes::calc")+
  paletteer::scale_fill_paletteer_d("ggthemes::calc")+
  My_Theme

ggsave(here::here("Output/Figures", "CALERIE_clocks_by_age_all_clocks.png"), width = 8, height = 8)

