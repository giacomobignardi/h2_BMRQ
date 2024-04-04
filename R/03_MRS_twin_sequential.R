#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Structural Equation Model Classical Twin Design specification : Trivariate sequential decomposition (Inspired by Cholesky)
#Program 03 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\XXX\04_packages')

# load packages
library(tidyverse)
library(lavaan)

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
twin_BMRQ_wide = twin_BMRQ_wide %>% rename(P1_T1= "music_discrimination_1", P2_T1 = "BAS_RR_1", P3_T1 = "BMRS_total_1", 
                          P1_T2= "music_discrimination_2", P2_T2 = "BAS_RR_2", P3_T2 = "BMRS_total_2") %>% 
  select(fam, zyg, age, P1_T1,P2_T1,P3_T1,P1_T2,P2_T2,P3_T2)

# SEQUENTIAL####
# source sequential model specification (Cholesky inspired but with Direct Symmetric Approach for the variance components)
source(sprintf("%s/%s/functions/lavaantwda/CTD_AE_sequential_5g_3p_lavaan.R", wdOA,wdOA_scripts))

Seq_AE_mod_fit = sem(Seq_AE_mod,
                          data = twin_BMRQ_wide,
                          group = "zyg",
                          group.label= c(1,2,3,4,5),
                          missing = "ML")

# Summary of sequential model
sumy_Seq_AE_mod_fit = summary(Seq_AE_mod_fit, fit.measures = TRUE, ci = T,standardized = TRUE)
fitMeasures(Seq_AE_mod_fit, c("CFI", "srmr"))

# final parameter estimates
sumy_Seq_AE_mod_fit$pe %>% filter(label %in% c("VarP1","VarP2","VarP3",
                                                    "h2_P1","h2_P2","h2_P3", "E_P3",
                                                    "varA31ut","varA32ut","varA33ut",
                                                    "a1", "a2", "a3",
                                                    "e1", "e2", "e3",
                                                    "R2a1","R2a2","R2a3","R2e1","R2e2","R2e3","R2ae23")) %>% 
  distinct(label, .keep_all = T)


# Supplementary file 6 ####
SFILE_6 = sumy_Seq_AE_mod_fit$pe %>% 
  distinct(label, .keep_all = T) %>% 
  select(label, est,std.all, se, pvalue, ci.lower,ci.upper) %>% 
  # remove defined parameters are the are obtained from the estimated paramteres
  slice(-c(30:nrow(sumy_Seq_AE_mod_fit$pe))) %>% 
  # remove path equal 1
  slice(-8) %>% 
  # remove covariances 
  filter(!(grepl("covMZ",label))) %>% 
  filter(!(grepl("covDZ",label))) %>% 
  mutate(label = recode(label,
                        "mean_f_P1" =  "M1w",
                        "mean_m_P1" =  "M1m",
                        "b1_P1" =  "B1",
                        "mean_f_P2" =  "M2w",
                        "mean_m_P2" =  "M2m",
                        "b1_P2" =  "B2",
                        "mean_f_P3" =  "M3w",
                        "mean_m_P3" =  "M3m",
                        "b1_P3" =  "B3",
                        "var_age" =  "σ2age",
                        "varA_P1" =  "σ2A1",
                        "varA_P2" =  "σ2A2",
                        "varA_P3" =  "σ2A3",
                        "varE_P1" =  "σ2E1",
                        "varE_P2" =  "σ2E2",
                        "varE_P3" =  "σ2E3",
                        "a1" =  "λA12",
                        "e1" =  "λE12",
                        "a2" =  "λA13",
                        "e2" =  "λE13",
                        "a3" =  "λA23",
                        "e3" =  "λE23",
                        
                        "mean_f_P3" =  "M3w",
                        "mean_m_P3" =  "M3m",
                        "b1_P3" =  "B3",
  )) %>% 
  rename(parameter = label,
         p = pvalue) %>% 
  mutate(across(where(is.numeric), round, 3))

write_csv(SFILE_6,sprintf("%s/%s/03_SFILE6.csv", wdOA,wdOA_output))

# P1 = music discrimination
# P2 = general reward sensitivity
# P3 = music reward sensitivity

# a1 = additive genetic path P1->P2
# a2 = additive genetic path P1->P3
# a3 = additive genetic path P2->P3

# e1 = residual path P1->P2
# e2 = residual path P1->P3
# e3 = residual path P2->P3

#  percentage of twin-h2 specific to music reward sensitivity
sumy_Seq_AE_mod_fit$pe %>% filter(label == "varAut") 

#  "deconfounded h2 music reward sensitivity
sumy_Seq_AE_mod_fit$pe %>% filter(label == "h2_P3") 

#  BMRQ variance explained
sumy_Seq_AE_mod_fit$pe %>% filter(grepl("R2a",label))
sumy_Seq_AE_mod_fit$pe %>% filter(grepl("R2e",label))

# VIZUALIZE####
## sequential decompostion of BMRQ variance####
sequential_data =  sumy_Seq_AE_mod_fit$pe%>% 
  filter(label %in% c("h2_P3", "E_P3", "R2a2", "R2e2", "R2a3", "R2e3")) %>% 
  distinct(label, .keep_all = T)

# save for external plotting
write_csv(sequential_data,sprintf("%s/%s/03_FIG2E_barplot_data.csv", wdOA,wdOA_output))

# proportion of BMRQ variance accounted for by each component
bar_plot = sumy_Seq_AE_mod_fit$pe%>% 
  filter(label %in% c("h2_P3", "E_P3", "R2a2", "R2e2", "R2a3", "R2e3")) %>% 
  distinct(label, .keep_all = T) %>%
  mutate(variable = "music reward sensitivity", 
         label = fct_recode(label,   
                            A3 = "h2_P3",
                            E3 = "E_P3",
                            A2= "R2a3",
                            E2= "R2e3",
                            A1 = "R2a2",
                            E1 = "R2e2")) %>% 
  mutate(component = factor(label, levels = rev(c("A1","E1","A2", "E2","A3","E3")))) %>% 
  ggplot(aes(variable,est, fill = component)) +
  geom_bar(position = "stack", stat = "identity", width = .33, color = "black") +
  scale_fill_manual(values = c("#5bbcd6","#f98400", "#00A08a","#F2ad00","#cccccc","#FD6467"))+
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.25),2), expand = c(0.0001,0.0001)) +
  geom_text(aes(y= c(.25,.4,.03,.08,NA,NA) , x = 1, label=round(est,2)), color="white", position = "stack") +
  labs(fill = "variance\ncomponent*",
       y = "proportion of variance") +
  theme_minimal(base_size = 12) +
theme(axis.ticks.x = element_blank())

# load phenotypic correlation
load(sprintf("%s/%s/02_pheno_cor_T1.Rdata",wdOA,wdOA_output))

#SAVE####
pdf(sprintf("%s/%s/03_Cholesky_estimates_BMRQ.pdf",wdOA,wdOA_ImageOutput),
    width =4.4,
    height =2.75)
bar_plot
dev.off()

#Fig 2cd
pdf(sprintf("%s/%s/03_Fig2cd_Cholesky.pdf",wdOA,wdOA_ImageOutput),
    width =8,
    height =5)
p_crosstrait_1|bar_plot
dev.off()


## Alternative estimator####
### MLR####
Seq_AE_mod_fit_MLR = sem(Seq_AE_mod,
                     data = twin_BMRQ_wide,
                     group = "zyg",
                     group.label= c(1,2,3,4,5),
                     missing = "ML",
                     estimator = "MLR")

# Summary of sequential model
sumy_Seq_AE_mod_fit_MLR = summary(Seq_AE_mod_fit_MLR, fit.measures = TRUE, ci = T)
fitMeasures(Seq_AE_mod_fit_MLR, c("CFI", "srmr"))

# standardized (over total variance including age) solution
sumy_std_Seq_AE_mod_fit = standardizedSolution(Seq_AE_mod_fit)
sumy_std_Seq_AE_mod_fit_MLR = standardizedSolution(Seq_AE_mod_fit_MLR)

# a path coefficients
sumy_Seq_AE_mod_fit_MLR$pe %>% filter(label == "a1" | label == "a2"| label == "a3") %>% distinct(label, .keep_all = T)
sumy_Seq_AE_mod_fit$pe %>% filter(label == "a1" | label == "a2"| label == "a3") %>% distinct(label, .keep_all = T)
sumy_std_Seq_AE_mod_fit %>%  filter(label == "a1" | label == "a2"| label == "a3") %>% distinct(label, .keep_all = T)
sumy_std_Seq_AE_mod_fit_MLR %>%  filter(label == "a1" | label == "a2"| label == "a3") %>% distinct(label, .keep_all = T)

# e path coefficients
sumy_Seq_AE_mod_fit_MLR$pe %>% filter(label == "e1" | label == "e2"| label == "e3") %>% distinct(label, .keep_all = T)
sumy_Seq_AE_mod_fit$pe %>% filter(label == "e1" | label == "e2"| label == "e3") %>% distinct(label, .keep_all = T)
sumy_std_Seq_AE_mod_fit %>%  filter(label == "e1" | label == "e2"| label == "e3") %>% distinct(label, .keep_all = T)
sumy_std_Seq_AE_mod_fit_MLR %>%  filter(label == "e1" | label == "e2"| label == "e3") %>% distinct(label, .keep_all = T)
