#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Strucutral Equation Model Classical Twin Design specification : Distinct factor solution vs Common genetic factor solution
#Program: 05 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\XXX\04_packages')

# load packages
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
twin_BMRQ_wide = twin_BMRQ_wide %>% rename(P1_T1= "BMRS_EE_1", P2_T1 = "BMRS_MR_1", P3_T1 = "BMRS_MS_1", P4_T1= "BMRS_SM_1", P5_T1 = "BMRS_SR_1",
                                           P1_T2= "BMRS_EE_2", P2_T2 = "BMRS_MR_2", P3_T2 = "BMRS_MS_2", P4_T2= "BMRS_SM_2", P5_T2 = "BMRS_SR_2") %>% 
  select(fam, zyg, age, P1_T1,P2_T1,P3_T1,P4_T1,P5_T1,P1_T2,P2_T2,P3_T2,P4_T2,P5_T2)

## SATURATED MULTIVARIATE MODEL####
# source saturated multivariate model specification
source(sprintf("%s/%s/functions/lavaantwda/CTD_SAT_5g_5p_lavaan.R", wdOA,wdOA_scripts))

# fit multivariate saturated model
SAT_fit_5v = sem(SAT_model_5v,
                          data = twin_BMRQ_wide,
                          group = "zyg",
                          group.label= c(1,2,3,4,5),
                          missing = "ML")

# summary stats of the saturated model
sumy_SAT_fit_5v = summary(SAT_fit_5v, fit.measures = TRUE, ci = T)

# fit indicies (note that the covariate age is manifests, so sat is partially constrained)
fitmeasures(SAT_fit_5v,c("cfi","srmr"))

# extract twin correlations stratified by zygosity
twin_cor_facets = sumy_SAT_fit_5v$pe %>% 
  filter(grepl("rp",label)) %>% #select correlations
  mutate(facet = substr(label,1,3),
         group = substr(label,5,nchar(label))) %>% 
  select(facet,group,est,ci.lower,ci.upper) %>%   
  mutate(group = fct_recode(group, MZf = "f_mz",
                            MZm = "m_mz",
                            DZf = "f_dz",
                            DZm = "m_dz",
                            DZos = "dos"),
         facet = fct_recode(facet, 
                            EE = "rp1",
                            MR = "rp2",
                            MS = "rp3",
                            SM = "rp4",
                            SR = "rp5")) %>% 
  mutate(group = factor(group, levels = rev(c("MZf","MZm","DZf","DZm","DZos")))) %>% 
  mutate(facet = factor(facet, levels = rev(c("EE" , "MR","MS" ,"SM", "SR"))))

## MULTIVARIATE AE MODEL####
source(sprintf("%s/%s/functions/lavaantwda/CTD_AE_CFShIPM_5g_5p.R", wdOA,wdOA_scripts))

## mDSA####
# fit correlated factor solution (multivarate via Direct Symmetric Approach)
correlated_AE_5g_fit = sem(correlated_AE_5g_model,
                 data = twin_BMRQ_wide,
                 group = "zyg",
                 group.label= c(1,2,3,4,5),
                 missing = "ML")

# save the model fit
save(correlated_AE_5g_fit,file = sprintf("%s/%s/correlated_AE_5g_fit.Rdata", wdOA,wdOA_output))

# summary of the mDSA
sumy_correlated_AE_5g_fit = summary(correlated_AE_5g_fit, ci = T, fit.measures = T, standardized = T)

#LRT: does the correlated factor solution worsen the saturated model fit?
lavTestLRT(SAT_fit_5v,correlated_AE_5g_fit)
fitmeasures(correlated_AE_5g_fit,c("cfi","srmr"))


# Supplementary file 7 ####
Sfile_7 = sumy_correlated_AE_5g_fit$pe%>% 
  distinct(label, .keep_all = T) %>% 
  select(label, est,std.all, se, pvalue, ci.lower,ci.upper) %>% 
  # remove defined parameters are the are obtained from the estimated paramteres
  slice(-c(78:nrow(sumy_correlated_AE_5g_fit$pe))) %>% 
  # remove path equal 1
  slice(-12) %>% 
  # remove covariances 
  filter(!(grepl("cov_MZ",label))) %>% 
  filter(!(grepl("cov_DZ",label))) %>% 
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
                        "mean_f_P4" =  "M4w",
                        "mean_m_P4" =  "M4m",
                        "b1_P4" =  "B4",
                        "mean_f_P5" =  "M5w",
                        "mean_m_P5" =  "M5m",
                        "b1_P5" =  "B5",
                        "var_age" =  "σ2age",
                        "varA_P1" =  "σ2A1",
                        "varA_P2" =  "σ2A2",
                        "varA_P3" =  "σ2A3",
                        "varA_P4" =  "σ2A4",
                        "varA_P5" =  "σ2A5",
                        "varE_P1" =  "σ2E1",
                        "varE_P2" =  "σ2E2",
                        "varE_P3" =  "σ2E3",
                        "varE_P4" =  "σ2E4",
                        "varE_P5" =  "σ2E5",
                        "covA_P12" =  "σA12",
                        "covA_P13" =  "σA13",
                        "covA_P14" =  "σA14",
                        "covA_P15" =  "σA15",
                        "covA_P23" =  "σA23",
                        "covA_P24" =  "σA24",
                        "covA_P25" =  "σA25",
                        "covA_P34" =  "σA34",
                        "covA_P35" =  "σA35",
                        "covA_P45" =  "σA45",
                        "covE_P12" =  "σE12",
                        "covE_P13" =  "σE13",
                        "covE_P14" =  "σE14",
                        "covE_P15" =  "σE15",
                        "covE_P23" =  "σE23",
                        "covE_P24" =  "σE24",
                        "covE_P25" =  "σE25",
                        "covE_P34" =  "σE34",
                        "covE_P35" =  "σE35",
                        "covE_P45" =  "σE45"
  )) %>% 
  rename(parameter = label,
         p = pvalue) %>% 
  mutate(across(where(is.numeric), round, 3))

write_csv(Sfile_7,sprintf("%s/%s/05_SFILE7.csv", wdOA,wdOA_output))


##hIPM####
# independent pathway model, hybrid common factor solution
ipm_hybrid_AE_5g_fit = sem(ipm_hybrid_AE_5g_model,
                           data = twin_BMRQ_wide,
                           group = "zyg",
                           group.label= c(1,2,3,4,5),
                           missing = "ML")

# summary of the hIPM
sumy_ipm_hybrid_AE_5g_fit = summary(ipm_hybrid_AE_5g_fit, ci = T, fit.measures = T, standardized = T)

# LRT: does the hybrid indepndent pathway model  worsen the saturated model fit?
lavTestLRT(correlated_AE_5g_fit ,ipm_hybrid_AE_5g_fit)
fitmeasures(ipm_hybrid_AE_5g_fit,c("cfi","srmr"))

# Supplementary file 7 ####
Sfile_8 = sumy_ipm_hybrid_AE_5g_fit$pe%>% 
  distinct(label, .keep_all = T) %>% 
  select(label, est,std.all, se, pvalue, ci.lower,ci.upper) %>% 
  # remove defined parameters are the are obtained from the estimated paramteres
  slice(-c(55:nrow(sumy_ipm_hybrid_AE_5g_fit$pe))) %>% 
  # remove path equal 1
  slice(-17) %>% 
  # remove covariances 
  filter(!(grepl("cov_MZ",label))) %>% 
  filter(!(grepl("cov_DZ",label))) %>% 
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
                        "mean_f_P4" =  "M4w",
                        "mean_m_P4" =  "M4m",
                        "b1_P4" =  "B4",
                        "mean_f_P5" =  "M5w",
                        "mean_m_P5" =  "M5m",
                        "b1_P5" =  "B5",
                        "var_age" =  "σ2age",
                        "A_P1" =  "λA1",
                        "A_P2" =  "λA2",
                        "A_P3" =  "λA3",
                        "A_P4" =  "λA4",
                        "A_P5" =  "λA5",
                        "varA_P1" =  "σ2Au1",
                        "varA_P2" =  "σ2Au2",
                        "varA_P3" =  "σ2Au3",
                        "varA_P4" =  "σ2Au4",
                        "varA_P5" =  "σ2Au5",
                        "varA_P1" =  "σ2Au1",
                        "varA_P2" =  "σ2Au2",
                        "varA_P3" =  "σ2Au3",
                        "varA_P4" =  "σ2Au4",
                        "varA_P5" =  "σ2Au5",
                        "varE_P1" =  "σ2E1",
                        "varE_P2" =  "σ2E2",
                        "varE_P3" =  "σ2E3",
                        "varE_P4" =  "σ2E4",
                        "varE_P5" =  "σ2E5",
                        "covE_P12" =  "σE12",
                        "covE_P13" =  "σE13",
                        "covE_P14" =  "σE14",
                        "covE_P15" =  "σE15",
                        "covE_P23" =  "σE23",
                        "covE_P24" =  "σE24",
                        "covE_P25" =  "σE25",
                        "covE_P34" =  "σE34",
                        "covE_P35" =  "σE35",
                        "covE_P45" =  "σE45"
  )) %>% 
  rename(parameter = label,
         p = pvalue) %>% 
  mutate(across(where(is.numeric), round, 3))

write_csv(Sfile_8,sprintf("%s/%s/05_SFILE8.csv", wdOA,wdOA_output))


## estimates####
# variance after accounting for age
sumy_correlated_AE_5g_fit$pe %>% filter(label == "varP1")

# standard deviation of the additive component
sumy_correlated_AE_5g_fit$pe %>% filter(label == "sdA_P1")

# heritability
sumy_correlated_AE_5g_fit$pe %>% filter(label == "h2_P1")

# biavariate heritability
sumy_correlated_AE_5g_fit$pe %>% filter(grepl("bh2", label) & grepl("1",label))

# genetic correlations
sumy_correlated_AE_5g_fit$pe%>% filter(grepl("rA", label) & !grepl("var",label)& grepl("1",label))

# residual correlations
sumy_correlated_AE_5g_fit$pe%>% filter(grepl("rE", label) & !grepl("var",label)& grepl("1",label))

# VIZUALIZE####
# correlation matrix####
# prapare the additive genetic lower diagonal
rA_matrix_row = sumy_correlated_AE_5g_fit$pe%>% filter(grepl("rA", label) & !grepl("var",label)) %>% 
  mutate(var1 = substr(lhs, nchar(lhs)-1,nchar(lhs)-1),
         var2 = substr(lhs, nchar(lhs),nchar(lhs))) %>% 
  select(var1,var2,est,ci.lower,ci.upper) %>% 
  mutate(est =round(est,2),
         ci.lower =round(ci.lower,2),
         ci.upper = round(ci.upper,2),
         var1 = recode(var1, "1" = "emotion\nevocation", 
                       "2" = "mood\nregulation", 
                       "3" = "music\nseeking", 
                       "4" = "sensory\nmotor", 
                       "5" = "social\nreward"),
         var2 = recode(var2, "1" = "emotion\nevocation", 
                       "2" = "mood\nregulation", 
                       "3" = "music\nseeking", 
                       "4" = "sensory\nmotor", 
                       "5" = "social\nreward"),
         type = "rA")

# prepare the unique environmental upper diagonal
rE_matrix_row = sumy_correlated_AE_5g_fit$pe%>% filter(grepl("rE", label) & !grepl("var",label)) %>% 
  mutate(var1 = substr(lhs, nchar(lhs)-1,nchar(lhs)-1),
         var2 = substr(lhs, nchar(lhs),nchar(lhs))) %>% 
  select(var1,var2,est,ci.lower,ci.upper) %>% 
  mutate(est =round(est,2),
         ci.lower =round(ci.lower,2),
         ci.upper = round(ci.upper,2),
         var1 = recode(var1, "1" = "emotion\nevocation", 
                       "2" = "mood\nregulation", 
                       "3" = "music\nseeking", 
                       "4" = "sensory\nmotor", 
                       "5" = "social\nreward"),
         var2 = recode(var2, "1" = "emotion\nevocation", 
                       "2" = "mood\nregulation", 
                       "3" = "music\nseeking", 
                       "4" = "sensory\nmotor", 
                       "5" = "social\nreward"),
         type = "rE")

# Make a copy of data but with the first two columns switched
rA_matrix_col = rA_matrix_row %>% select(2:1, 3:6)
rE_matrix_col = rE_matrix_row %>% select(2:1, 3:6)
names(rA_matrix_col) = names(rA_matrix_row)
names(rE_matrix_col) = names(rE_matrix_row)

# Stick the two data frames together
rA_matrix = rbind(rA_matrix_row, rA_matrix_col)
rE_matrix = rbind(rE_matrix_row, rE_matrix_col)

# romove redundand information
rA_matrix <- rA_matrix[which(as.numeric(as.factor(rA_matrix$var2)) >
                               as.numeric(as.factor(rA_matrix$var1))),]

rE_matrix <- rE_matrix[which(as.numeric(as.factor(rE_matrix$var2)) <
                               as.numeric(as.factor(rE_matrix$var1))),]

# Stick the two data frames together.
r_matrix <- rbind(rA_matrix, rE_matrix)

# Create the confidence intervals using paste
r_matrix$CI <- paste0("(", r_matrix$ci.lower, ", ", r_matrix$ci.upper, ")")

# select heritability
h_diagonal = sumy_correlated_AE_5g_fit$pe %>% 
  filter(grepl("h2",label) & !grepl("bh2",label)) %>%
  mutate(var1 = seq(1:5),var2 = seq(1:5))%>% 
  select(var1,var2,est,ci.lower,ci.upper) %>% 
  mutate(est =round(est,2),
         ci.lower =round(ci.lower,2),
         ci.upper = round(ci.upper,2),
         var1 = recode(var1, "1" = "emotion\nevocation", 
                       "2" = "mood\nregulation", 
                       "3" = "music\nseeking", 
                       "4" = "sensory\nmotor", 
                       "5" = "social\nreward"),
         var2 = recode(var2, "1" = "emotion\nevocation", 
                       "2" = "mood\nregulation", 
                       "3" = "music\nseeking", 
                       "4" = "sensory\nmotor", 
                       "5" = "social\nreward"),
         type = "h2")

# Create the heritability text by using paste
h_diagonal$CI <- paste0("(", h_diagonal$ci.lower, ", ", h_diagonal$ci.upper, ")")

# Stick the two data frames together.
r_h_matrix <- rbind(r_matrix ,h_diagonal)

# convert to factor so that rows and columns have the same order as the data
r_h_matrix = r_h_matrix %>% 
  mutate(var1 = factor(var1, levels = c("emotion\nevocation","mood\nregulation","music\nseeking","sensory\nmotor","social\nreward" )),
         var2 = factor(var2, levels = rev(c("emotion\nevocation","mood\nregulation","music\nseeking","sensory\nmotor","social\nreward"))))

# save for external plotting
write_csv(r_h_matrix,sprintf("%s/%s/05_FIG3E_matrix_data.csv", wdOA,wdOA_output))

# final summary matrix of the correlated factor solution results
cor_plot = 
  ggplot(mapping = aes(var1,var2)) +  
  coord_fixed()+
  geom_tile(data = r_h_matrix %>% filter(type == "rE"), aes(fill = est), color = "white")+
  geom_text(data = r_h_matrix %>% filter(type == "rE"), aes(label = est), size = 5, position = position_nudge(y = 0.1),color="white") +
  geom_text(data = r_h_matrix %>% filter(type == "rE"), aes(label = CI), size = 3.5, position = position_nudge(y = -0.1),color="white") +
  scale_fill_gradient2(na.value = "white", high = "#022e39",
                       low = "#044855", limits = c(0,1))+
  labs(fill = "rE")+
  new_scale_fill() +
  geom_tile(data = r_h_matrix %>% filter(type == "rA"), aes(fill = est), color = "white")+
  geom_text(data = r_h_matrix %>% filter(type == "rA"), aes(label = est), size = 5, position = position_nudge(y = 0.1),color="white") +
  geom_text(data = r_h_matrix %>% filter(type == "rA"), aes(label = CI), size = 3.5, position = position_nudge(y = -0.1),color="white") +
  scale_fill_gradient2(na.value = "white", high = "#c0424e",
                       low = "#d94b58", limits = c(0,1))+
  labs(fill = "rA")+
  new_scale_fill() +
  geom_tile(data = r_h_matrix %>% filter(type == "h2"), aes(fill = est), color = "white")+
  geom_text(data = r_h_matrix %>% filter(type == "h2"), aes(label = est), size = 5, position = position_nudge(y = 0.1),color="white") +
  geom_text(data = r_h_matrix %>% filter(type == "h2"), aes(label = CI), size = 3.5, position = position_nudge(y = -0.1),color="white") +
  scale_fill_gradient2(na.value = "white", high = "#793f4e",
                       low = "white", limits = c(0,1))+
  labs(x = "",
       y = "",
       fill = expression(h^2)) +
  theme_classic(base_size = 12)+
  scale_y_discrete(limits=rev)+
  theme(axis.text       = element_text(size = 12),
        axis.title      = element_text(size = 12), 
        text            = element_text(size = 12),
        axis.line = element_blank(),
        legend.position = "bottom")

#noe on the average pair-wise genetic correlation per facet 
rbind(
  ee = r_h_matrix %>% filter(type == "rA" & var1 == "emotion\nevocation") %>% mutate(est_z = psych::fisherz(est)) %>% summarise(mean =  psych::fisherz2r(mean(est_z)), sd = psych::fisherz2r(sd(est_z))),
  mr = r_h_matrix %>% filter((type == "rA" & var1 == "mood\nregulation")  | (type == "rA" & var2 == "mood\nregulation")) %>% mutate(est_z = psych::fisherz(est)) %>% summarise(mean =  psych::fisherz2r(mean(est_z)), sd = psych::fisherz2r(sd(est_z))),
  ms = r_h_matrix %>% filter((type == "rA" & var1 == "music\nseeking")  | (type == "rA" & var2 == "music\nseeking")) %>% mutate(est_z = psych::fisherz(est)) %>% summarise(mean =  psych::fisherz2r(mean(est_z)), sd = psych::fisherz2r(sd(est_z))),
  sm = r_h_matrix %>% filter((type == "rA" & var1 == "sensory\nmotor")  | (type == "rA" & var2 == "sensory\nmotor")) %>% mutate(est_z = psych::fisherz(est)) %>% summarise(mean =  psych::fisherz2r(mean(est_z)), sd = psych::fisherz2r(sd(est_z))),
  sr = r_h_matrix %>% filter((type == "rA" & var1 == "social\nreward")  | (type == "rA" & var2 == "social\nreward")) %>% mutate(est_z = psych::fisherz(est)) %>% summarise(mean =  psych::fisherz2r(mean(est_z)), sd = psych::fisherz2r(sd(est_z)))
)

#SAVE####
pdf(sprintf("%s/%s/05_rArEh2_cor_plot.pdf",wdOA,wdOA_ImageOutput),
    height = 7,
    width =7)
cor_plot
dev.off()

## Alternative estimator####
### MLR####
# mDSA
# fit correlated factor solution (multivarate via Direct Symmetric Approach)
correlated_AE_5g_fit_MLR = sem(correlated_AE_5g_model,
                           data = twin_BMRQ_wide,
                           group = "zyg",
                           group.label= c(1,2,3,4,5),
                           missing = "ML",
                           estimator = "MLR")
# summary of the mDSA
sumy_correlated_AE_5g_fit = summary(correlated_AE_5g_fit_MLR, ci = T, fit.measures = T)

# hIPM#
# independent pathway model, hybrid common factor solution
ipm_hybrid_AE_5g_fit_MLR = sem(ipm_hybrid_AE_5g_model,
                           data = twin_BMRQ_wide,
                           group = "zyg",
                           group.label= c(1,2,3,4,5),
                           missing = "ML",
                           estimator = "MLR")

# summary of the hIPM
sumy_ipm_hybrid_AE_5g_fit_MLR = summary(ipm_hybrid_AE_5g_fit_MLR, ci = T, fit.measures = T)

# LRT: does the hybrid indepndent pathway model  worsen the saturated model fit?
lavTestLRT(correlated_AE_5g_fit ,ipm_hybrid_AE_5g_fit)
lavTestLRT(correlated_AE_5g_fit_MLR ,ipm_hybrid_AE_5g_fit_MLR)