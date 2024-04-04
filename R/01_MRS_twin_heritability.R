#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Strucutral Equation Modeling Classical Twin Design specification
#Program: 01------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\XXX\04_packages')

# load packages
library(tidyverse)
library(lavaan)
library(patchwork)

# clean working environment 
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
source(sprintf("%s/%s/functions/lavaantwda/CTD_SAT_5g_univariate_lavaan.R", wdOA,wdOA_scripts))

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
qsumy_SAT_fit = summary(qSAT_fit,ci = T, fit.measures = T)

# fit indicies (note that the covariate age is manifests, so sat is partially constrained)
fitmeasures(qSAT_fit,c("cfi","srmr"))

# extract standardized covariances (cor) for twin pairs
sumy_SATstd_fit = standardizedsolution(qSAT_fit)
sumy_SATstd_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_f_mz")
sumy_SATstd_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_m_mz")
sumy_SATstd_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_f_dz")
sumy_SATstd_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_m_dz")
sumy_SATstd_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_dos")

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
rbind(
  lavTestLRT(qSAT_fit,SAT_fit_bor_2),
  lavTestLRT(qSAT_fit,SAT_fit_zyg_3),
  lavTestLRT(qSAT_fit,SAT_fit_sxm_4),
  lavTestLRT(qSAT_fit,SAT_fit_sxv_5),
  lavTestLRT(qSAT_fit,SAT_fit_qnt_6),
  lavTestLRT(qSAT_fit,SAT_fit_qal_7)
)

# summary of model after testing for assumptions
sumy_SATstd_final_fit = standardizedsolution(SAT_fit_qal_7)
sumy_SAT_final_fit = summary(SAT_fit_qal_7, ci = T)
summary(qSAT_fit)

# means
sumy_SAT_final_fit$pe %>% filter(label == "mean_f")
sumy_SAT_final_fit$pe %>% filter(label == "mean_m")

# standardized b1
sumy_SATstd_final_fit %>% filter(label == "b1")

# which variance component do we expect? if rMZ is twice higher than rDZ then ADE, otherwise ACE
rMZ = sumy_SATstd_final_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_mz") %>% slice(1)
rDZ = sumy_SATstd_final_fit %>% select(label,est.std,ci.lower,ci.upper) %>% filter(label == "cov_dz") %>% slice(1)
expected  = ifelse(rMZ$est.std > 2*rDZ$est.std, "ADE", "ACE")
# We expect A(D)E model

# ADE MODEL####
# #ADE model####
# source ADE model specification
source(sprintf("%s/%s/functions/lavaantwda/CTD_ADE_5g_univariate_lavaan.R", wdOA,wdOA_scripts))

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
sumy_ADE_fit = summary(ADE_fit,ci = T, fit.measures = T)
sumy_ADE_fit

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
fitMeasures(AE_fit,c("cfi","srmr"))

# summary statistics final model
sumy_AE_fit = summary(AE_fit,ci = T, fit.measures = T)

# tidy results ADE
ADE_ADE_fit = sumy_ADE_fit$pe %>%
  filter(label == "A" | label == "D" | label == "E") %>% 
  select(label, est) %>% 
  mutate(model = "ADE",
         component = factor(label, levels = c("E", "D", "A")))

# tidy results AE
ADE_AE_fit = sumy_AE_fit$pe %>% 
  filter(label == "A" | label == "E") %>% 
  add_row(label = "D") %>% select(label, est) %>% 
  mutate(model = "AE",
         component = factor(label, levels = c("E", "D", "A")))

# Supplementary file 5 ####
# prepare table for publication
SFile_5 = sumy_AE_fit$pe %>% 
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

write_csv(SFile_5,sprintf("%s/%s/01_SFILE5.csv", wdOA,wdOA_output))

# VIZUALIZE####
## phenotypic correlations####
# extract correlations and confidence intervals from saturated model with covariate age
pheno_cor = sumy_SATstd_fit %>% select(label,est.std,ci.lower,ci.upper) %>% 
  mutate(label = fct_recode(label,   
                              MZf = "cov_f_mz",
                              MZm = "cov_m_mz",
                              DZf = "cov_f_dz",
                              DZm = "cov_m_dz",
                              DZos= "cov_dos")) %>% 
  mutate(zygosity = factor(label, levels = rev(c("MZf","MZm","DZf","DZm","DZos"))))

# save for external plotting
write_csv(pheno_cor,sprintf("%s/%s/01_FIG1B_forestplot_data.csv", wdOA,wdOA_output))

# plot phenotypic correlations and confidence intervals stratified by zygosity
forest_like_plot = pheno_cor %>% 
  filter(zygosity %in% c("MZf", "MZm", "DZf", "DZm", "DZos")) %>% 
  ggplot(aes(est.std,zygosity)) +
  geom_point(position = position_dodge(.9), size = 2.5)+#3182BD
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper), width = .1, size = 1,position = position_dodge(.9)) +
  xlim(-0.01,1) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  scale_color_viridis_d()+
  labs(
    y = "zygosity",
    x = "phenotypic r"
    )+
  theme_classic(12)+
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")

## A(D)E model estimates####
ADE_estimates = rbind(ADE_ADE_fit,ADE_AE_fit)

# save for external plotting
write_csv(ADE_estimates,sprintf("%s/%s/01_FIG1C_barplot_data.csv", wdOA,wdOA_output))

# plot ADE and AE model estimates
bar_plot = rbind(ADE_ADE_fit,ADE_AE_fit) %>% 
  mutate(model = factor(model, levels = c("ADE", "AE"))) %>%
  rename(estimate = est) %>% 
  ggplot(aes(model,estimate, fill = component)) +
  geom_bar(stat = "identity", width = .5, color = "black") +
  scale_fill_manual(values = c("#5bbcd6","#B19D61" ,"#f98400" ))+
  geom_text(aes(x=c(1,1,1,2,2,2), y = c(.21,.28,.3,.27,0,.50), label=round(estimate,2)), 
            color="white", position = "stack")+
  labs(x = "model",
       y="proportion of variance")+
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.25),2), expand = c(0.001,0.001)) +
  theme_minimal()+
  theme(axis.ticks.x = element_blank())

#SAVE####
pdf(sprintf("%s/%s/01_twin_correlations_BMRQ.pdf",wdOA,wdOA_ImageOutput),
    width =3.5,
    height =2.75)
forest_like_plot
dev.off()

pdf(sprintf("%s/%s/01_ADE_estimates_BMRQ.pdf",wdOA,wdOA_ImageOutput),
    width =4.4,
    height =2.75)
bar_plot
dev.off()

#F1bc
pdf(sprintf("%s/%s/01_Fig1bc_ADE_BMRQ.pdf",wdOA,wdOA_ImageOutput),
    width =4,
    height =5)
forest_like_plot/bar_plot
dev.off()

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

# model check (info on test are in the CTD_SAT_5g_lavaan.R script)
SAT_fit_sxm_4_MLR = sem(SAT_mod_sxm_4, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML", estimator ="MLR")

# model comparison: sex means
lavTestLRT(qSAT_fit,SAT_fit_sxm_4)
lavTestLRT(qSAT_fit_MLR,SAT_fit_sxm_4_MLR)

# parameter: age
SAT_fit_qal_7_MLR = sem(SAT_mod_qal_7, data = twin_BMRQ_wide, group = "zyg", group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")
sumy_SATstd_final_fit_MLR = standardizedsolution(SAT_fit_qal_7_MLR)

# standardized b1
sumy_SATstd_final_fit_MLR %>% filter(label == "b1")
sumy_SATstd_final_fit %>% filter(label == "b1")


# AE model
# fit parsimonious models by dropping D and E
AE_fit_MLR = sem(AE_model,data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")
E_fit_MLR = sem(E_model,data = twin_BMRQ_wide,group = "zyg",group.label= c(1,2,3,4,5), missing = "ML", estimator = "MLR")

## model comparison####
lavTestLRT(AE_fit_MLR,E_fit_MLR)

fitMeasures(ADE_fit,c("cfi","srmr"))
fitMeasures(AE_fit,c("cfi","srmr"))


