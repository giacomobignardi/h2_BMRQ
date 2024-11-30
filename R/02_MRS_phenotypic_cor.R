#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Compute observed phenotypic correlations between music discrimination abilities, general reward sensitivity, and music reward sensitivity
#Program: 02 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\C4_Behavior_Genectics_Unit\\Giacomo Bignardi\\2023_BMRQ\\04_packages')

# load packages
library(tidyverse)

# clean working environment 
rm(list = ls())
set.seed(123)

# set Open Access working directories
wdOA = getwd()
wdOA_scripts = "06_github/R"
wdOA_output = "03_outputs/01_final"
wdOA_sfile = "06_github/SFILE"

# set not Open Access working directories
wdNOA_Data = "01_data"
wdOA_ImageOutput = "05_images/01_final"

# load dataFrames:
twin_BMRQ_wide  = read_csv(sprintf("%s/%s/00_twin_BMRQ_wide.csv", wdOA,wdOA_output))

# CROSS-TRAIT CORRELATIONS####
# twin 1####
phenotypes_t1 = twin_BMRQ_wide %>% select(music_discrimination_1,BAS_RR_1,BMRS_total_1) %>% rename("music\nperceptual\nabilities" = music_discrimination_1, "general\nreward\nsensitivty" = BAS_RR_1,"music\nreward\nsensitivty" = BMRS_total_1)
observed_varCov_T1 =  cor(phenotypes_t1, use = "pairwise.complete.obs")

# compute standard deviations (Sd)
sd_mda_t1 = sd(phenotypes_t1$`music\nperceptual\nabilities`, na.rm = T)
sd_grs_t1 = sd(phenotypes_t1$`general\nreward\nsensitivty`, na.rm = T)
sd_mrs_t1 = sd(phenotypes_t1$`music\nreward\nsensitivty`, na.rm = T)

# substitute sd to diagonal values
diag(observed_varCov_T1) = c(sd_mda_t1,sd_grs_t1,sd_mrs_t1)

# prepare dataframe
observed_varCov_T1[lower.tri(observed_varCov_T1, diag = F)] = NA 
observed_varCov_T1 = observed_varCov_T1 %>% 
  as.data.frame() %>% 
  rownames_to_column()

# adapted from StackOverlfow https://stackoverflow.com/questions/67530405/change-orientation-of-diagonal-of-correlation-plot-using-ggcorrplot-package-if
observed_varCov_long_T1 = observed_varCov_T1 %>% 
  pivot_longer(c("music\nperceptual\nabilities": "music\nreward\nsensitivty"),names_to = "colname", values_to = "value" )

# convert to factor so that rows and columns have the same order as the data
observed_varCov_long_T1 = observed_varCov_long_T1 %>% 
  mutate(rowname = factor(rowname, levels = c("music\nperceptual\nabilities", "general\nreward\nsensitivty","music\nreward\nsensitivty")),
         colname = factor(colname, levels = rev(c("music\nperceptual\nabilities", "general\nreward\nsensitivty","music\nreward\nsensitivty"))))

# set diagonal and the top-right half of the matrix to 0 so that those cells appears white
p_crosstrait_1 = ggplot(observed_varCov_long_T1, aes(x = colname, y = rowname, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "white", high = "#3B4371",
                       low = "#F3904F",mid = "#EAE1DF",
                       limits = c(-0.5,0.5)) +
  geom_text(aes(colname, rowname, label = round(value,2)), color = "black", size = 3.5) +
  theme_minimal(base_size = 12)+
  labs(x = NULL, y = NULL, fill =  "cross\nphenotypic r")+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_equal(1)

# twin 2####
phenotypes_T2 = twin_BMRQ_wide %>% select(music_discrimination_2,BAS_RR_2,BMRS_total_2) %>% rename("music\nperceptual\nabilities" = music_discrimination_2, "general\nreward\nsensitivty" = BAS_RR_2,"music\nreward\nsensitivty" = BMRS_total_2)
observed_varCov_T2 =  cor(phenotypes_T2, use = "pairwise.complete.obs")

# compute standard deviations (Sd)
sd_mda_T2 = sd(phenotypes_T2$`music\nperceptual\nabilities`, na.rm = T)
sd_grs_T2 = sd(phenotypes_T2$`general\nreward\nsensitivty`, na.rm = T)
sd_mrs_T2 = sd(phenotypes_T2$`music\nreward\nsensitivty`, na.rm = T)

# substitute sd to diagonal values
diag(observed_varCov_T2) = c(sd_mda_T2,sd_grs_T2,sd_mrs_T2)

observed_varCov_T2[lower.tri(observed_varCov_T2, diag = F)] = NA 
observed_varCov_T2 = observed_varCov_T2 %>% 
  as.data.frame() %>% 
  rownames_to_column()

# adapted from StackOverlfow https://stackoverflow.com/questions/67530405/change-orientation-of-diagonal-of-correlation-plot-using-ggcorrplot-package-if
observed_varCov_long_T2 = observed_varCov_T2 %>% 
  pivot_longer(c("music\nperceptual\nabilities": "music\nreward\nsensitivty"),names_to = "colname", values_to = "value" )

# convert to factor so that rows and columns have the same order as the data
observed_varCov_long_T2 = observed_varCov_long_T2 %>% 
  mutate(rowname = factor(rowname, levels = c("music\nperceptual\nabilities", "general\nreward\nsensitivty","music\nreward\nsensitivty")),
         colname = factor(colname, levels = rev(c("music\nperceptual\nabilities", "general\nreward\nsensitivty","music\nreward\nsensitivty"))))

# set diagonal and the top-right half of the matrix to 0 so that those cells appears white
p_crosstrait_2 = ggplot(observed_varCov_long_T2, aes(x = colname, y = rowname, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "white", high = "#3B4371",
                       low = "#F3904F",mid = "#EAE1DF",
                       limits = c(-0.5,0.5)) +
  geom_text(aes(colname, rowname, label = round(value,2)), color = "black", size = 3.5) +
  theme_minimal()+
  labs(x = NULL, y = NULL, fill = "cross\nphenotypic r")+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_equal(1)

#SAVE####
pdf(sprintf("%s/%s/02_pheno_cor.pdf",wdOA,wdOA_ImageOutput),
    width =7,
    height =3.5)
(p_crosstrait_1 + theme(legend.position = "none")|p_crosstrait_2) + plot_layout(guides ="collect") + plot_annotation(tag_levels = "a")
dev.off()