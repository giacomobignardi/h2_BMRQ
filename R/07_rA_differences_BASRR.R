#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Strucutral Equation Model Classical Twin Design specification : Distinct Factor Solution including general reward
#Program: 07 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\XXX\04_packages')

# load packages
.libPaths('W:\\C4_Behavior_Genectics_Unit\\Giacomo Bignardi\\2023_BMRQ\\04_packages')
library(tidyverse)
library(lavaan)
library(wesanderson)
library(patchwork)
library(ggnewscale)

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
twin_BMRQ_wide = twin_BMRQ_wide %>% rename(P1_T1= "BAS_RR_1", P2_T1 = "BMRS_EE_1", P3_T1 = "BMRS_MR_1", P4_T1 = "BMRS_MS_1", P5_T1= "BMRS_SM_1", P6_T1 = "BMRS_SR_1",
                                           P1_T2= "BAS_RR_2", P2_T2 = "BMRS_EE_2", P3_T2 = "BMRS_MR_2", P4_T2 = "BMRS_MS_2", P5_T2= "BMRS_SM_2", P6_T2 = "BMRS_SR_2") %>% 
  select(fam, zyg, age, P1_T1,P2_T1,P3_T1,P4_T1,P5_T1,P6_T1,P1_T2,P2_T2,P3_T2,P4_T2,P5_T2,P6_T2)

## MULTIVARIATE AE MODEL####
source(sprintf("%s/%s/functions/lavantwda/CTD_AE_CFShIPM_5g_6p.R", wdOA,wdOA_scripts))

## mDSA####
# fit correlated factor solution (multivarate via Direct Symmetric Approach)
correlated_AE_6g_fit = sem(correlated_AE_6g_model,
                 data = twin_BMRQ_wide,
                 group = "zyg",
                 group.label= c(1,2,3,4,5),
                 missing = "ML")

# summary of the mDSA
sumy_correlated_AE_6g_fit = summary(correlated_AE_6g_fit, ci = T)

# fit statisitcs of the mDSA
fitmeasures(correlated_AE_6g_fit, c("cfi","srmr"))

# rA
rA__correlated_AE_6g_fit =  sumy_correlated_AE_6g_fit$pe %>% filter(grepl("rA",lhs) & grepl("rA_P1",label))
rE__correlated_AE_6g_fit =  sumy_correlated_AE_6g_fit$pe %>% filter(grepl("rE",lhs) & grepl("rA_P1",label))

# save rA for external plotting
write_csv(rA__correlated_AE_6g_fit,sprintf("%s/%s/07_rA_BASRR.csv", wdOA,wdOA_output))

# differences between covA
sumy_correlated_AE_6g_fit$pe %>% filter(grepl("d_",lhs) & pvalue<.05)

sumy_correlated_AE_6g_fit$pe %>% filter(grepl("covA",label) )

# sanity check
ggplot(rA__correlated_AE_6g_fit, aes(est,lhs))+
  geom_point() +
  geom_pointrange(aes(xmin = ci.lower, xmax = ci.upper))+
  theme_minimal()


