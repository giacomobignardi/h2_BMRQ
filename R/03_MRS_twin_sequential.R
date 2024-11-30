#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Structural Equation Model Classical Twin Design specification : Trivariate sequential decomposition (Inspired by Cholesky)
#Program 03 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\C4_Behavior_Genectics_Unit\\Giacomo Bignardi\\2023_BMRQ\\04_packages')

# load packages
library(tidyverse)
library(lavaan)

# clean working enviroment 
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
twin_BMRQ_wide = twin_BMRQ_wide %>% rename(P1_T1= "music_discrimination_1", P2_T1 = "BAS_RR_1", P3_T1 = "BMRS_total_1", 
                                           P1_T2= "music_discrimination_2", P2_T2 = "BAS_RR_2", P3_T2 = "BMRS_total_2") %>% 
  select(fam, zyg, age, P1_T1,P2_T1,P3_T1,P1_T2,P2_T2,P3_T2)

# SEQUENTIAL####
# source sequential model specification (Cholesky inspired but with Direct Symmetric Approach for the variance components)
source(sprintf("%s/%s/functions/lavaantwda/AE_sequential_5g_3p.R", wdOA,wdOA_scripts))

Seq_AE_fit = sem(Seq_AE_mod,
                 data = twin_BMRQ_wide,
                 group = "zyg",
                 group.label= c(1,2,3,4,5),
                 missing = "ML")

# Summary of sequential model
Seq_AE_sumy = summary(Seq_AE_fit, fit.measures = TRUE, ci = T,standardized = TRUE)
# Summary of Cholesky
Seq_AE_std = standardizedsolution(Seq_AE_fit)
fitMeasures(Seq_AE_fit, c("CFI", "srmr"))

# final parameter estimates
Seq_AE_sumy$pe %>% filter(label %in% c("VarP1","VarP2","VarP3",
                                       "h2_P1","h2_P2","h2_P3", 
                                       "e2_P1","e2_P2","e2_P3",
                                       "varA31ut","varA32ut","varA33ut",
                                       "a12", "a13", "a23",
                                       "e12", "e13", "e23",
                                       "R2a12","R2a13","R2a23","R2e12","R2e13","R2e23","R2ae123")) %>% 
  
  distinct(label, .keep_all = T)

# LRT: significance testing
Seq_AE_fit_noa12 = sem(paste(Seq_AE_mod,lrt_a12),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
Seq_AE_fit_noa13 = sem(paste(Seq_AE_mod,lrt_a13),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
Seq_AE_fit_noa23 = sem(paste(Seq_AE_mod,lrt_a23),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
Seq_AE_fit_noe12 = sem(paste(Seq_AE_mod,lrt_e12),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
Seq_AE_fit_noe13 = sem(paste(Seq_AE_mod,lrt_e13),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
Seq_AE_fit_noe23 = sem(paste(Seq_AE_mod,lrt_e23),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")

lavTestLRT(Seq_AE_fit ,Seq_AE_fit_noa12)
lavTestLRT(Seq_AE_fit ,Seq_AE_fit_noa13)
lavTestLRT(Seq_AE_fit ,Seq_AE_fit_noa23)
lavTestLRT(Seq_AE_fit ,Seq_AE_fit_noe12)
lavTestLRT(Seq_AE_fit ,Seq_AE_fit_noe13)
lavTestLRT(Seq_AE_fit ,Seq_AE_fit_noe23)

# Supplementary File 2 ####
SFILE_2 = Seq_AE_std %>% 
  mutate(label = ifelse(lhs == "A_P1_T1" & rhs == "P1_T1", "a1",label),
         label = ifelse(lhs == "A_P2_T1" & rhs == "P2_T1", "a2",label),
         label = ifelse(lhs == "A_P3_T1" & rhs == "P3_T1", "a3",label),
         label = ifelse(lhs == "E_P1_T1" & rhs == "P1_T1", "e1",label),
         label = ifelse(lhs == "E_P2_T1" & rhs == "P2_T1", "e2",label),
         label = ifelse(lhs == "E_P3_T1" & rhs == "P3_T1", "e3",label)) %>% 
  distinct(label, .keep_all = T) %>% 
  select(label, est.std, se, pvalue, ci.lower,ci.upper)  %>% 
  filter(!(grepl("covMZ",label))) %>% 
  filter(!(grepl("covDZ",label))) %>% 
  filter(!(grepl("h2",label))) %>% 
  filter(!(grepl("e1_",label))) %>% 
  filter(!(grepl("e2_",label))) %>% 
  filter(!(grepl("e3_",label))) %>% 
  filter(!(grepl("var",label))) %>% 
  filter(est.std != 0) %>% 
  filter(label != "") %>%
  # remove defined parameters are the are obtained from the estimated paramteres
  slice(1:21) %>% 
  # remove path equal 1
  # remove covariances 
  filter(!(grepl("covMZ",label))) %>% 
  filter(!(grepl("covDZ",label))) %>% 
  mutate(label = recode(label,
                        "mean_f_P1" =  "M1w",
                        "mean_m_P1" =  "M1m",
                        "b1_P1" =  "β1",
                        "mean_f_P2" =  "M2w",
                        "mean_m_P2" =  "M2m",
                        "b1_P2" =  "β2",
                        "mean_f_P3" =  "M3w",
                        "mean_m_P3" =  "M3m",
                        "b1_P3" =  "β3",
                        "var_age" =  "σ2age",
                        "a1" =  "λA1",
                        "e1" =  "λE1",
                        "a2" =  "λA2",
                        "e2" =  "λE2",
                        "a3" =  "λA3",
                        "e3" =  "λE3",
                        "a12" =  "λA12",
                        "e12" =  "λE12",
                        "a13" =  "λA13",
                        "e13" =  "λE13",
                        "a23" =  "λA23",
                        "e23" =  "λE23"
  )) %>% 
  rename(parameter = label,
         p = pvalue) %>% 
  mutate(rad_est.sdt = ifelse(startsWith(parameter,"λ"),est.std^2,NA))%>% 
  mutate(across(where(is.numeric), round, 3))%>% 
  arrange(parameter)

write_csv(SFILE_2,sprintf("%s/%s/03_supplementary_data_2.csv", wdOA,wdOA_sfile))

# P1 = music discrimination
# P2 = general reward sensitivity
# P3 = music reward sensitivity

# a12 = additive genetic path P1->P2
# a13 = additive genetic path P1->P3
# a23 = additive genetic path P2->P3

# e12 = residual path P1->P2
# e13 = residual path P1->P3
# e23 = residual path P2->P3

#  percentage of twin-h2 specific to music reward sensitivity
Seq_AE_sumy$pe %>% filter(label == "varA33ut") 

#  adjusted h2 music reward sensitivity
Seq_AE_sumy$pe %>% filter(label == "h2_P3") 

#  BMRQ variance explained
Seq_AE_sumy$pe %>% filter(grepl("R2a",label))
Seq_AE_sumy$pe %>% filter(grepl("R2e",label))

# VIZUALIZE####
## sequential decompostion of BMRQ variance####
sequential_data =  Seq_AE_sumy$pe%>% 
  filter(label %in% c("h2_P3", "E_P3", "R2a2", "R2e2", "R2a3", "R2e3")) %>% 
  distinct(label, .keep_all = T)

# save for external plotting
write_csv(sequential_data,sprintf("%s/%s/03_FIG2_barplot_data.csv", wdOA,wdOA_output))

## Alternative estimator####
### MLR####
Seq_AE_fit_MLR = sem(Seq_AE_mod,
                     data = twin_BMRQ_wide,
                     group = "zyg",
                     group.label= c(1,2,3,4,5),
                     missing = "ML",
                     estimator = "MLR")

# Summary of sequential model
Seq_AE_sumy_MLR = summary(Seq_AE_fit_MLR, fit.measures = TRUE, ci = T,standardized = TRUE)

# Summary of Cholesky
Seq_AE_std_MLR = standardizedsolution(Seq_AE_fit_MLR)

# to report in Supplementary table 5 (pMLR)
Seq_AE_std_MLR %>% 
  distinct(label, .keep_all = T) %>% 
  select(label, est.std, se, pvalue, ci.lower,ci.upper)  %>% 
  filter(!(grepl("covMZ",label))) %>% 
  filter(!(grepl("covDZ",label))) %>% 
  filter(!(grepl("var",label))) %>% 
  filter(est.std != 0) %>% 
  # remove defined parameters are the are obtained from the estimated paramteres
  slice(1:21) %>% 
  # remove path equal 1
  # remove covariances 
  filter(!(grepl("covMZ",label))) %>% 
  filter(!(grepl("covDZ",label))) %>% 
  mutate(label = recode(label,
                        "mean_f_P1" =  "M1w",
                        "mean_m_P1" =  "M1m",
                        "b1_P1" =  "β1",
                        "mean_f_P2" =  "M2w",
                        "mean_m_P2" =  "M2m",
                        "b1_P2" =  "β2",
                        "mean_f_P3" =  "M3w",
                        "mean_m_P3" =  "M3m",
                        "b1_P3" =  "β3",
                        "var_age" =  "σ2age",
                        "a1" =  "λA1",
                        "e1" =  "λE1",
                        "a2" =  "λA2",
                        "e2" =  "λE2",
                        "a3" =  "λA3",
                        "e3" =  "λE3",
                        "a12" =  "λA12",
                        "e12" =  "λE12",
                        "a13" =  "λA13",
                        "e13" =  "λE13",
                        "a23" =  "λA23",
                        "e23" =  "λE23"
  )) %>% 
  rename(parameter = label,
         p = pvalue) %>% 
  mutate(rad_est.sdt = ifelse(startsWith(parameter,"λ"),est.std^2,NA))%>% 
  mutate(across(where(is.numeric), round, 3))%>% 
  arrange(parameter)
