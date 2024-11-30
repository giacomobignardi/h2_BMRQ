#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Strucutral Equation Modeling Classical Twin Design specification
#Program: 01------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\C4_Behavior_Genectics_Unit\\Giacomo Bignardi\\2023_BMRQ\\04_packages')

# load packages
library(tidyverse)
library(lavaan)
library(patchwork)

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

# Supplemental scatter-plot of correlation 
scatter_plot = twin_BMRQ_wide %>% 
  mutate(zyg = recode_factor(
    zyg,
    `1` = "MZw",
    `2` = "MZm",
    `3` = "DZw",
    `4` = "DZm",
    `5` = "DZos",
  )) %>% 
  ggplot(aes(BMRS_total_1,BMRS_total_2)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_grid(cols = vars(zyg)) +
  labs(x = "music reward sensitivity T1",
       y = "music reward sensitivity T2",
       color = "zygosity")+
  theme_bw()

#SAVE####
pdf(sprintf("%s/%s/01_scatter_plot_BMRQ.pdf",wdOA,wdOA_ImageOutput),
    width =9,
    height =2.5)
scatter_plot
dev.off()

# SATURATED MODEL####
# source SAT model specification
source(sprintf("%s/%s/functions/lavaantwda/SAT_5g_univariate.R", wdOA,wdOA_scripts))

# rename to meet lavaan custom code requirements
twin_BMRQ_wide = twin_BMRQ_wide %>% 
  rename(P_T1 = BMRS_total_1,
         P_T2 = BMRS_total_2)


# sanity check that there are no duplicated families
n_distinct(twin_BMRQ_wide$fam)
nrow(twin_BMRQ_wide)

# fit saturated model with age as covariate
qSAT_fit = sem(qSAT_mod,
               data = twin_BMRQ_wide,
               group = "zyg",
               group.label= c(1,2,3,4,5),
               missing = "ML"
)

# sanity check: group order is respected
qSAT_fit@Data@group.label

# summary saturated model
qSAT_sumy = summary(qSAT_fit,ci = T, fit.measures = T)

# fit indicies (note that the covariate age is manifests, so sat is partially constrained)
fitmeasures(qSAT_fit,c("cfi","srmr"))
#cfi  srmr 
#1.000 0.008 

# extract standardized covariances (cor) for twin pairs
qSAT_std = standardizedsolution(qSAT_fit)
qSAT_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_f_mz")
qSAT_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_m_mz")
qSAT_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_f_dz")
qSAT_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_m_dz")
qSAT_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_dos")

# model check (info on test are in the CTD_SAT_5g_lavaan.R script)
SAT_fit = sem(SAT_mod, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_cov_0 = sem(SAT_mod_cov_0, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_cov_1 = sem(SAT_mod_cov_1, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_bor_2 = sem(SAT_mod_bor_2, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_zyg_3 = sem(SAT_mod_zyg_3, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_sxm_4 = sem(SAT_mod_sxm_4, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_sxv_5 = sem(SAT_mod_sxv_5, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_qnt_6 = sem(SAT_mod_qnt_6, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")
SAT_fit_qal_7 = sem(SAT_mod_qal_7, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML")

## model comparison####
# saturated vs baseline: reminder, here we constrain age var across groups to ease interpretations of standardized variances across groups
lavTestLRT(SAT_fit,qSAT_fit)

rbind(
  lavTestLRT(qSAT_fit,SAT_fit_bor_2),
  lavTestLRT(qSAT_fit,SAT_fit_zyg_3),
  lavTestLRT(qSAT_fit,SAT_fit_sxm_4),
  lavTestLRT(qSAT_fit,SAT_fit_sxv_5),
  lavTestLRT(qSAT_fit,SAT_fit_qnt_6),
  lavTestLRT(qSAT_fit,SAT_fit_qal_7)
)

# report test against saturated (included already in lavaan output "Model Test User Model")
SAT_fit_bor_2
SAT_fit_zyg_3
SAT_fit_sxm_4
SAT_fit_sxv_5
SAT_fit_qnt_6
SAT_fit_qal_7

# test for covariate age (LRT addition)
## equality of paths
lavTestLRT(SAT_fit_cov_0,SAT_fit_cov_1)

# summary of model after testing for assumptions
SAT_final_std = standardizedsolution(SAT_fit_qal_7)
SAT_final_sumy = summary(SAT_fit_qal_7, ci = T)
summary(qSAT_fit)

# means
SAT_final_sumy$pe %>% filter(label == "mean_f")
SAT_final_sumy$pe %>% filter(label == "mean_m")

# Coehn s d
SAT_final_sumy$pe %>% filter(label == "d")
# alternative d (uncostrained variances)
SAT_zyg_3_sumy = summary(SAT_fit_zyg_3)
SAT_zyg_3_sumy$pe %>% filter(label == "d") # same to the rounding

# standardized b1
SAT_final_std %>% filter(label == "b1")
# unstandardized b1
SAT_final_sumy$pe %>% filter(label == "b1")

# which variance component do we expect? if rMZ is twice higher than rDZ then ADE, otherwise ACE
rMZ = SAT_final_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_mz") %>% slice(1)
rDZ = SAT_final_std %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_dz") %>% slice(1)
expected  = ifelse(rMZ$est.std > 2*rDZ$est.std, "ADE", "ACE")
# We expect A(D)E model

# ADE MODEL####
# #ADE model####
# source ADE model specification
source(sprintf("%s/%s/functions/lavaantwda/ADE_5g_univariate.R", wdOA,wdOA_scripts))

# fit ADE model with age as covariate
ADE_fit = sem(ADE_model,
              data = twin_BMRQ_wide,
              group = "zyg",
              group.label= c(1,2,3,4,5),
              missing = "ML"
)

# sanity check; group order is respected
SAT_fit@Data@group.label

# summary ADE model
ADE_sumy = summary(ADE_fit,ci = T, fit.measures = T)
ADE_sumy

## AE model####
# fit parsimonious models by dropping D and E
AE_fit = sem(AE_model,data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5), missing = "ML")
E_fit = sem(E_model,data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5), missing = "ML")

## model comparison####
rbind(
  lavTestLRT(SAT_fit,ADE_fit),
  lavTestLRT(ADE_fit,AE_fit),
  lavTestLRT(AE_fit,E_fit)
)

fitMeasures(ADE_fit,c("cfi","srmr"))
#0.983 0.041

fitMeasures(AE_fit,c("cfi","srmr"))
#0.981 0.042

# summary statistics final model
AE_sumy = summary(AE_fit,ci = T, fit.measures = T)

# tidy results ADE
ADE = ADE_sumy$pe %>%
  filter(label == "A" | label == "D" | label == "E") %>% 
  select(label, est) %>% 
  mutate(model = "ADE",
         component = factor(label, levels = c("E", "D", "A")))

# tidy results AE
AE = AE_sumy$pe %>% 
  filter(label == "A" | label == "E") %>% 
  add_row(label = "D") %>% select(label, est) %>% 
  mutate(model = "AE",
         component = factor(label, levels = c("E", "D", "A")))

# Supplementary file 1 ####
# prepare table for publication
SFile_1 = AE_sumy$pe %>% 
  distinct(label, .keep_all = T) %>% 
  select(label, est, se, pvalue, ci.lower,ci.upper) %>% 
  # remove path equal 1
  slice(-4) %>% 
  # remove covariances (as they are equal to the alpha*sigma2A)
  filter(!(label  == "covMZA" | label  == "covDZA" )) %>% 
  # remove defined parameters are the are obtained from the estimated paramteres
  slice(-c(7:12)) %>% 
  mutate(label = recode(label,
                        "mean_f" =  "Mw",
                        "mean_m" =  "Mm",
                        "b1" =  "Bage",
                        "var_age" =  "σ2age",
                        "varA" =  "σ2A",
                        "varE" =  "σ2E",
  )) %>% 
  rename(parameter = label,
         p = pvalue) %>% 
  mutate(across(where(is.numeric), round, 3))

write_csv(SFile_1,sprintf("%s/%s/01_supplementary_data_1.csv", wdOA,wdOA_sfile))

# VIZUALIZE####
## phenotypic correlations####
qSAT_std = standardizedsolution(qSAT_fit)
# extract correlations and confidence intervals from saturated model with covariate age
pheno_cor = qSAT_std  %>% select(label,est.std,ci.lower,ci.upper) %>% 
  mutate(label = fct_recode(label,   
                            MZf = "cov_f_mz",
                            MZm = "cov_m_mz",
                            DZf = "cov_f_dz",
                            DZm = "cov_m_dz",
                            DZos= "cov_dos")) %>% 
  mutate(zygosity = factor(label, levels = rev(c("MZf","MZm","DZf","DZm","DZos"))))

# save for external plotting
write_csv(pheno_cor,sprintf("%s/%s/01_FIG1_forestplot_data.csv", wdOA,wdOA_output))

## A(D)E model estimates####
ADE_estimates = rbind(ADE,AE)

# save for external plotting
write_csv(ADE_estimates,sprintf("%s/%s/01_FIG1_barplot_data.csv", wdOA,wdOA_output))

## Alternative estimator####
### MLR####
# fit saturated model with age as covariate
qSAT_fit_MLR = sem(qSAT_mod,
                   data = twin_BMRQ_wide,
                   group = "zyg",
                   group.label= c(1,2,3,4,5),
                   missing = "ML",
                   estimator = "MLR",
)


# fit indicies (note that the covariate age is manifests, so sat is partially constrained)
fitmeasures(qSAT_fit_MLR,c("cfi","srmr"))
#cfi  srmr 
#1.000 0.008 

# model check (info on test are in the CTD_SAT_5g_lavaan.R script)
SAT_fit_sxm_4_MLR = sem(SAT_mod_sxm_4, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML", estimator ="MLR")

# model comparison: sex means
lavTestLRT(qSAT_fit,SAT_fit_sxm_4)
lavTestLRT(qSAT_fit_MLR,SAT_fit_sxm_4_MLR)

# parameter: age
SAT_fit_cov_0_MLR = sem(SAT_fit_cov_0, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")
SAT_fit_cov_1_MLR = sem(SAT_mod_cov_1, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")

## equality of paths
lavTestLRT(SAT_fit_cov_0_MLR,SAT_fit_cov_1_MLR)

# AE model
# fit parsimonious models by dropping D and E
AE_fit_MLR = sem(AE_model,data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")
E_fit_MLR = sem(E_model,data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")

## model comparison####
lavTestLRT(AE_fit_MLR,E_fit_MLR)