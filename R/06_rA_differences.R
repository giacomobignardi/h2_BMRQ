#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Strucutral Equation Model Classical Twin Design specification : Distinct Factor Solution including music perceptual ability
#Program: 06 ------------------------------------------------------------------------------------------------------------------------------
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
twin_BMRQ_wide = twin_BMRQ_wide %>% rename(P1_T1= "music_discrimination_1", P2_T1 = "BMRS_EE_1", P3_T1 = "BMRS_MR_1", P4_T1 = "BMRS_MS_1", P5_T1= "BMRS_SM_1", P6_T1 = "BMRS_SR_1",
                                           P1_T2= "music_discrimination_2", P2_T2 = "BMRS_EE_2", P3_T2 = "BMRS_MR_2", P4_T2 = "BMRS_MS_2", P5_T2= "BMRS_SM_2", P6_T2 = "BMRS_SR_2") %>% 
  select(fam, zyg, age, P1_T1,P2_T1,P3_T1,P4_T1,P5_T1,P6_T1,P1_T2,P2_T2,P3_T2,P4_T2,P5_T2,P6_T2)

## MULTIVARIATE AE MODEL####
source(sprintf("%s/%s/functions/lavaantwda/AE_hIPM_5g_6p.R", wdOA,wdOA_scripts))

## mDSA####
# fit correlated factor solution (multivarate via Direct Symmetric Approach)
cfs_AE_fit = sem(cfs_AE_mod,
                 data = twin_BMRQ_wide,
                 group = "zyg",
                 group.label= c(1,2,3,4,5),
                 missing = "ML")

# summary of the mDSA
cfs_AE_sumy = summary(cfs_AE_fit, ci = T)

# fit statisitcs of the mDSA
fitmeasures(cfs_AE_fit, c("cfi","srmr"))

# rA
rA_AE =  cfs_AE_sumy$pe %>% filter(grepl("rA",lhs) & grepl("rA_P1",label) )
rE_AE =  cfs_AE_sumy$pe %>% filter(grepl("rE",lhs) & grepl("rE_P1",label) )
cA_AE =  cfs_AE_sumy$pe %>% filter(grepl("covA",label)) %>% distinct(label,.keep_all = T)

# save rA for external plotting
write_csv(rA_AE,sprintf("%s/%s/06_rA_SMDT.csv", wdOA,wdOA_output))
write_csv(rE_AE,sprintf("%s/%s/06_rE_SMDT.csv", wdOA,wdOA_output))

# differences between back transfomed zA
cfs_AE_sumy$pe %>% 
  filter(grepl("dzA_",lhs)) %>% 
  mutate(rest = psych::fisherz2r(est))

cfs_AE_sumy$pe %>% 
  filter(grepl("dzE_",lhs)) %>% 
  mutate(rest = psych::fisherz2r(est))
sig_rA = cfs_AE_sumy$pe %>% filter(grepl("dzA_",lhs) & pvalue<.05 ) %>% pull(lhs)
cfs_AE_sumy$pe %>% filter(grepl("rA",label))

# sanity check
ggplot(rA_AE, aes(est,lhs))+
  geom_point() +
  geom_pointrange(aes(xmin = ci.lower, xmax = ci.upper))+
  theme_minimal()

# LRT based signficance test
# Compute restricited models
AP123_sig_fit = sem(paste(cfs_AE_mod,constrain_AP123),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP124_sig_fit = sem(paste(cfs_AE_mod,constrain_AP124),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP125_sig_fit = sem(paste(cfs_AE_mod,constrain_AP125),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP126_sig_fit = sem(paste(cfs_AE_mod,constrain_AP126),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP134_sig_fit = sem(paste(cfs_AE_mod,constrain_AP134),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP135_sig_fit = sem(paste(cfs_AE_mod,constrain_AP135),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP136_sig_fit = sem(paste(cfs_AE_mod,constrain_AP136),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP145_sig_fit = sem(paste(cfs_AE_mod,constrain_AP145),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP146_sig_fit = sem(paste(cfs_AE_mod,constrain_AP146),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
AP156_sig_fit = sem(paste(cfs_AE_mod,constrain_AP156),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")

# apply model comparison
lavTestLRT(cfs_AE_fit,AP123_sig_fit)
lavTestLRT(cfs_AE_fit,AP124_sig_fit)
lavTestLRT(cfs_AE_fit,AP125_sig_fit)
lavTestLRT(cfs_AE_fit,AP126_sig_fit)
lavTestLRT(cfs_AE_fit,AP134_sig_fit)
lavTestLRT(cfs_AE_fit,AP135_sig_fit)
lavTestLRT(cfs_AE_fit,AP136_sig_fit)
lavTestLRT(cfs_AE_fit,AP145_sig_fit)
lavTestLRT(cfs_AE_fit,AP146_sig_fit)
lavTestLRT(cfs_AE_fit,AP156_sig_fit)

# Compute restricited models
EP123_sig_fit = sem(paste(cfs_AE_mod,constrain_EP123),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP124_sig_fit = sem(paste(cfs_AE_mod,constrain_EP124),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP125_sig_fit = sem(paste(cfs_AE_mod,constrain_EP125),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP126_sig_fit = sem(paste(cfs_AE_mod,constrain_EP126),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP134_sig_fit = sem(paste(cfs_AE_mod,constrain_EP134),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP135_sig_fit = sem(paste(cfs_AE_mod,constrain_EP135),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP136_sig_fit = sem(paste(cfs_AE_mod,constrain_EP136),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP145_sig_fit = sem(paste(cfs_AE_mod,constrain_EP145),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP146_sig_fit = sem(paste(cfs_AE_mod,constrain_EP123),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")
EP156_sig_fit = sem(paste(cfs_AE_mod,constrain_EP156),data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5),missing = "ML")

# apply model comparison
lavTestLRT(cfs_AE_fit,EP123_sig_fit)
lavTestLRT(cfs_AE_fit,EP124_sig_fit)
lavTestLRT(cfs_AE_fit,EP125_sig_fit)
lavTestLRT(cfs_AE_fit,EP126_sig_fit)
lavTestLRT(cfs_AE_fit,EP134_sig_fit)
lavTestLRT(cfs_AE_fit,EP135_sig_fit)
lavTestLRT(cfs_AE_fit,EP136_sig_fit)
lavTestLRT(cfs_AE_fit,EP145_sig_fit)
lavTestLRT(cfs_AE_fit,EP146_sig_fit)
lavTestLRT(cfs_AE_fit,EP156_sig_fit)
