#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Compute observed phenotypic correlations between BMRQ facet
#Program: 04 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\XXX\04_packages')

# load packages
library(tidyverse)

# clean working enviroment 
rm(list = ls())
set.seed(123)

# set Open Access working directories
wdOA = getwd()
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs"

# set not Open Access working directories
wdNOA_Data = "01_data"
wdOA_ImageOutput = "05_images"

# load dataFrames:
twin_BMRQ_wide  = read_csv(sprintf("%s/%s/00_twin_BMRQ_wide.csv", wdOA,wdOA_output))

# CROSS-TRAIT CORRELATIONS####
# twin 1####
observed_varCov_T1 =  cor(twin_BMRQ_wide  %>% 
                            rename("emotion\nevocation" = BMRS_EE_1, "mood\nregulation" = BMRS_MR_1,"music\nseeking" = BMRS_MS_1, "sensory\nmotor" = BMRS_SM_1,"social\nreward" = BMRS_SR_1) %>% 
                            select("emotion\nevocation","mood\nregulation","music\nseeking","sensory\nmotor","social\nreward" ), use = "pairwise.complete.obs")

# weaker pairwise correlation
cor.test(twin_BMRQ_wide$BMRS_MS_1,twin_BMRQ_wide$BMRS_SM_1)
# stronger pairwise correlation
cor.test(twin_BMRQ_wide$BMRS_EE_1,twin_BMRQ_wide$BMRS_MR_1)

# compute standard deviations (Sd)
sd_ee_t1 = sd(twin_BMRQ_wide$BMRS_EE_1, na.rm = T)
sd_mr_t1 = sd(twin_BMRQ_wide$BMRS_MR_1, na.rm = T)
sd_ms_t1 = sd(twin_BMRQ_wide$BMRS_MS_1, na.rm = T)
sd_sm_t1 = sd(twin_BMRQ_wide$BMRS_SM_1, na.rm = T)
sd_sr_t1 = sd(twin_BMRQ_wide$BMRS_SR_1, na.rm = T)

# substitute sd to diagonal values
diag(observed_varCov_T1) = c(sd_ee_t1,sd_mr_t1,sd_ms_t1,sd_sm_t1,sd_sr_t1)

observed_varCov_T1[upper.tri(observed_varCov_T1, diag = F)] = NA 
observed_varCov_T1 = observed_varCov_T1 %>% 
  as.data.frame() %>% 
  rownames_to_column()

# adapted from StackOverlfow https://stackoverflow.com/questions/67530405/change-orientation-of-diagonal-of-correlation-plot-using-ggcorrplot-package-if
observed_varCov_long_T1 = observed_varCov_T1 %>% 
  pivot_longer(c("emotion\nevocation":"social\nreward"),names_to = "colname", values_to = "value" )

# convert to factor so that rows and columns have the same order as the data
observed_varCov_long_T1 = observed_varCov_long_T1 %>% 
  mutate(rowname = factor(rowname, levels = rev(c(unique(rowname)))),
         colname = factor(colname, levels = c(unique(colname))))

# set diagonal and the top-right half of the matrix to 0 so that those cells appears white
p_crossfacet_1 = ggplot(observed_varCov_long_T1, aes(x = colname, y = rowname, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "white", high = "#3B4371",
                       low = "#F3904F", limits = c(0,1)) +
  geom_text(aes(colname, rowname, label = round(value,2)), color = "black", size = 3.5) +
  theme_minimal()+
  labs(x = NULL, y = NULL, subtitle = "",fill = "phenotypic\ncross-facet r")+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_equal(1)

# twin 2####
observed_varCov_T2 =  cor(twin_BMRQ_wide  %>% 
                            rename("emotion\nevocation" = BMRS_EE_2, "mood\nregulation" = BMRS_MR_2,"music\nseeking" = BMRS_MS_2, "sensory\nmotor" = BMRS_SM_2,"social\nreward" = BMRS_SR_2) %>% 
                            select("emotion\nevocation","mood\nregulation","music\nseeking","sensory\nmotor","social\nreward" ), use = "pairwise.complete.obs")

# compute standard deviations (Sd)
sd_ee_T2 = sd(twin_BMRQ_wide$BMRS_EE_2, na.rm = T)
sd_mr_T2 = sd(twin_BMRQ_wide$BMRS_MR_2, na.rm = T)
sd_ms_T2 = sd(twin_BMRQ_wide$BMRS_MS_2, na.rm = T)
sd_sm_T2 = sd(twin_BMRQ_wide$BMRS_SM_2, na.rm = T)
sd_sr_T2 = sd(twin_BMRQ_wide$BMRS_SR_2, na.rm = T)

# substitute sd to diagonal values
diag(observed_varCov_T2) = c(sd_ee_T2,sd_mr_T2,sd_ms_T2,sd_sm_T2,sd_sr_T2)

observed_varCov_T2[upper.tri(observed_varCov_T2, diag = F)] = NA 
observed_varCov_T2 = observed_varCov_T2 %>% 
  as.data.frame() %>% 
  rownames_to_column()

# adapted from StackOverlfow https://stackoverflow.com/questions/67530405/change-orientation-of-diagonal-of-correlation-plot-using-ggcorrplot-package-if
observed_varCov_long_T2 = observed_varCov_T2 %>% 
  pivot_longer(c("emotion\nevocation":"social\nreward"),names_to = "colname", values_to = "value" )

# convert to factor so that rows and columns have the same order as the data
observed_varCov_long_T2 = observed_varCov_long_T2 %>% 
  mutate(rowname = factor(rowname, levels = rev(c(unique(rowname)))),
         colname = factor(colname, levels = c(unique(colname))))

# set diagonal and the top-right half of the matrix to 0 so that those cells appears white
p_crossfacet_2 = ggplot(observed_varCov_long_T2, aes(x = colname, y = rowname, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(na.value = "white", high = "#3B4371",
                       low = "#F3904F", limits = c(0,1)) +
  geom_text(aes(colname, rowname, label = round(value,2)), color = "black", size = 3.5) +
  theme_minimal()+
  labs(x = NULL, y = NULL, subtitle = "",fill = "phenotypic\ncross-facet r")+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_equal(1)


#SAVE####
pdf(sprintf("%s/%s/04_pheno_cor_facet_T1.pdf",wdOA,wdOA_ImageOutput),
    width =4,
    height =3.5)
p_crossfacet_1
dev.off()

pdf(sprintf("%s/%s/04_pheno_cor_facet_T2.pdf",wdOA,wdOA_ImageOutput),
    width =4,
    height =3.5)
p_crossfacet_2
dev.off()