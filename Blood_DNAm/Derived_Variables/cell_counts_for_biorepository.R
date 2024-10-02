pacman::p_load(tidyverse, patchwork)

# Cell counts for Biorepository

cells <- read_csv("../Data_Resources_FOR_BIOREPOSITORY/Derived_Variables/CALERIE_CPR_Blood_12_EPIC_Extended_Cell_Counts.csv")
clocks <-read_csv("../Data_Resources_FOR_BIOREPOSITORY/Derived_Variables//CALERIE_agedwb_horvath_datout.output_cleaned_extra_variables.csv") %>% select(barcode, Female, Age, CR, fu)

df <-left_join(clocks, cells, by = 'barcode') %>% filter(fu == 'base') %>% select(barcode:Treg)

all_long <-pivot_longer(df, cols = -c(barcode), names_to = "Cell Type", values_to = "Cell Proportion")

cols <-   paletteer::paletteer_d("ggthemes::Classic_Cyclic") 


#####################################################################################
#
#####################################################################################
  
  # Extract non_neutrophil data
  non_neu_long <-all_long %>% 
    filter(`Cell Type` != "Neu")
  # And take factor levels
  my_names <-levels(as_factor(non_neu_long$`Cell Type`))
  
  # Extract neutrophil data
  neu_long <-all_long %>% 
    filter(`Cell Type` == "Neu")
  
  # Plot the data
  # Plot non-neutrophil data
  non_neu_plot <-non_neu_long %>% 
    # Remove na for categorical variable
    ggplot() +
    stat_boxplot(aes(x = `Cell Type`, y = `Cell Proportion`, 
                     color = `Cell Type`),
                 geom = "errorbar",
                 lwd = 1) +
    geom_jitter(width = 0.1, alpha = 0.4, aes(x = `Cell Type`, y = `Cell Proportion`, fill = `Cell Type`), pch = 21, color = 'gray50') +
    geom_boxplot(aes(x = `Cell Type`, y = `Cell Proportion`, 
                     color = `Cell Type`),
                 outlier.shape = NA,
                 lwd = 1, alpha = 0.3) +
    theme_bw(base_size = 14) + 
    labs(x = "Cell Types", y = "DNAm-estimated Proportion of White Blood Cells") + 
    theme(legend.position = "top") +
    scale_x_discrete(labels = c(my_names)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5, size = 14)) +
    theme(axis.title.x = element_blank())+
    scale_fill_manual(values = c(cols))+
    scale_color_manual(values = c(cols))
  
  # Plot neutrophil data
  neu_plot <-neu_long %>% 
    ggplot() +
    stat_boxplot(aes(x = `Cell Type`, y = `Cell Proportion`, 
                     color = `Cell Type`),
                 geom = "errorbar",
                 lwd = 1) +
    geom_jitter(width = 0.1, alpha = 0.4, aes(x = `Cell Type`, y = `Cell Proportion`), pch = 21, color = 'gray50', fill = cols[12]) +
    geom_boxplot(aes(x = `Cell Type`, y = `Cell Proportion`, 
                    color = `Cell Type`),
                 outlier.shape = NA,
                 lwd = 1, alpha = 0.3) +
    theme_bw(base_size = 14) + 
    labs(x = "Cell Types", y = "DNAm-estimated Proportion of White Blood Cells") + 
    theme(legend.position = "top") +
    scale_x_discrete(labels = "Neutrophils") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5, size = 14)) +
    theme(axis.title.x = element_blank())+
    scale_color_manual(values = cols[12])

  
  
    # Join and print the figure itself
non_neu_plot + theme(legend.position = "none")+
    neu_plot + theme(legend.position = "none") + plot_layout(ncol = 2, widths = c(9.3, 1)) 

ggsave(here::here("Output/Figures", "Cell_Boxplots_for_Biorepository.png"), width = 8, height = 6)











