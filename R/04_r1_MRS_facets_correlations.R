#Author: Giacomo Bignardi
#
#
#
#Description: obtain overarching summary of pairwise correlations stratified across zygosity groups in response to Reviewer 1 Revision 1
#Program: 00_descriptive_summary ------------------------------------------------------------------------------------------------------------------------------
# load packages
.libPaths('W:\\C4_Behavior_Genectics_Unit\\Giacomo Bignardi\\2023_BMRQ\\04_packages')
library(tidyverse)
library(openxlsx)

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
twin_BMRQ_wide_var = twin_BMRQ_wide %>% select(zyg,
                          music_discrimination_1,music_discrimination_2,
                          BAS_RR_1,     BAS_RR_2,
                          BMRS_total_1, BMRS_total_2,
                          BMRS_EE_1,    BMRS_EE_2,
                          BMRS_MR_1,    BMRS_MR_2,
                          BMRS_MS_1,    BMRS_MS_2,
                          BMRS_SM_1,    BMRS_SM_2,
                          BMRS_SR_1,    BMRS_SR_2) 

# compute cross trait cross twin correlations
cor_MZw = cor(twin_BMRQ_wide_var %>% filter(zyg ==1) %>% select(-c(zyg)), use = "pairwise.complete.obs")
cor_MZm = cor(twin_BMRQ_wide_var %>% filter(zyg ==2) %>% select(-c(zyg)), use = "pairwise.complete.obs")
cor_DZw = cor(twin_BMRQ_wide_var %>% filter(zyg ==3) %>% select(-c(zyg)), use = "pairwise.complete.obs")
cor_DZm = cor(twin_BMRQ_wide_var %>% filter(zyg ==4) %>% select(-c(zyg)), use = "pairwise.complete.obs")
cor_DZos = cor(twin_BMRQ_wide_var %>% filter(zyg ==5) %>% select(-c(zyg)), use = "pairwise.complete.obs")

# compute complete sample for cross trait cross twin correlations
np_MZw = psych::pairwiseCount(twin_BMRQ_wide_var %>% filter(zyg ==1) %>% select(-c(zyg)))
np_MZm = psych::pairwiseCount(twin_BMRQ_wide_var %>% filter(zyg ==2) %>% select(-c(zyg)))
np_DZw = psych::pairwiseCount(twin_BMRQ_wide_var %>% filter(zyg ==3) %>% select(-c(zyg)))
np_DZm = psych::pairwiseCount(twin_BMRQ_wide_var %>% filter(zyg ==4) %>% select(-c(zyg)))
np_DZos =psych::pairwiseCount(twin_BMRQ_wide_var %>% filter(zyg ==5) %>% select(-c(zyg)))

# Reorder rows and columns to have  phenotypic correlations in the upper left and within- and cross- trait correlations in the bottom left
row_order_cor_MZw = order(grepl("2$", rownames(cor_MZw)), grepl("1$", rownames(cor_MZw)))
col_order_cor_MZw = order(grepl("2$", colnames(cor_MZw)), grepl("1$", colnames(cor_MZw)))
cor_MZw = cor_MZw[row_order_cor_MZw, row_order_cor_MZw]

row_order_cor_MZm = order(grepl("2$", rownames(cor_MZm)), grepl("1$", rownames(cor_MZm)))
col_order_cor_MZm = order(grepl("2$", colnames(cor_MZm)), grepl("1$", colnames(cor_MZm)))
cor_MZm = cor_MZm[row_order_cor_MZm, row_order_cor_MZm]

row_order_cor_DZw = order(grepl("2$", rownames(cor_DZw)), grepl("1$", rownames(cor_DZw)))
col_order_cor_DZw = order(grepl("2$", colnames(cor_DZw)), grepl("1$", colnames(cor_DZw)))
cor_DZw = cor_DZw[row_order_cor_DZw, row_order_cor_DZw]

row_order_cor_DZm = order(grepl("2$", rownames(cor_DZm)), grepl("1$", rownames(cor_DZm)))
col_order_cor_DZm = order(grepl("2$", colnames(cor_DZm)), grepl("1$", colnames(cor_DZm)))
cor_DZm = cor_DZm[row_order_cor_DZm, row_order_cor_DZm]

row_order_cor_DZos = order(grepl("2$", rownames(cor_DZos)), grepl("1$", rownames(cor_DZos)))
col_order_cor_DZos = order(grepl("2$", colnames(cor_DZos)), grepl("1$", colnames(cor_DZos)))
cor_DZos = cor_DZos[row_order_cor_DZos, row_order_cor_DZos]


# Reorder rows and columns to have  phenotypic npairs in the upper left and within- and cross- trait correlations in the bottom left
row_order_np_MZw = order(grepl("2$", rownames(np_MZw)), grepl("1$", rownames(np_MZw)))
col_order_np_MZw = order(grepl("2$", colnames(np_MZw)), grepl("1$", colnames(np_MZw)))
np_MZw = np_MZw[row_order_np_MZw, row_order_np_MZw]

row_order_np_MZm = order(grepl("2$", rownames(np_MZm)), grepl("1$", rownames(np_MZm)))
col_order_np_MZm = order(grepl("2$", colnames(np_MZm)), grepl("1$", colnames(np_MZm)))
np_MZm = np_MZm[row_order_np_MZm, row_order_np_MZm]

row_order_np_DZw = order(grepl("2$", rownames(np_DZw)), grepl("1$", rownames(np_DZw)))
col_order_np_DZw = order(grepl("2$", colnames(np_DZw)), grepl("1$", colnames(np_DZw)))
np_DZw = np_DZw[row_order_np_DZw, row_order_np_DZw]

row_order_np_DZm = order(grepl("2$", rownames(np_DZm)), grepl("1$", rownames(np_DZm)))
col_order_np_DZm = order(grepl("2$", colnames(np_DZm)), grepl("1$", colnames(np_DZm)))
np_DZm = np_DZm[row_order_np_DZm, row_order_np_DZm]

row_order_np_DZos = order(grepl("2$", rownames(np_DZos)), grepl("1$", rownames(np_DZos)))
col_order_np_DZos = order(grepl("2$", colnames(np_DZos)), grepl("1$", colnames(np_DZos)))
np_DZos = np_DZos[row_order_np_DZos, row_order_np_DZos]


# rename 
phen_names =c(
  "SMDT_tot_T1",
  "BAS_rr_T1",
  "BMRQ_tot_T1",
  "emo_evo_f1_T1",
  "mod_reg_f2_T1",
  "mus_sek_f3_T1",
  "sen_mtr_f4_T1",
  "soc_rwd_f5_T1",
  "SMDT_tot_T2",
  "BAS_rr_T2",
  "BMRQ_tot_T2",
  "emo_evo_f1_T2",
  "mod_reg_f2_T2",
  "mus_sek_f3_T2",
  "sen_mtr_f4_T2",
  "soc_rwd_f5_T2"
)
colnames(cor_MZw) = phen_names; rownames(cor_MZw) = phen_names 
colnames(cor_MZm) = phen_names; rownames(cor_MZm) = phen_names
colnames(cor_DZw) = phen_names; rownames(cor_DZw) = phen_names 
colnames(cor_DZm) = phen_names; rownames(cor_DZm) = phen_names
colnames(cor_DZos)= phen_names; rownames(cor_DZos) = phen_names

colnames(np_MZw) = phen_names; rownames(np_MZw) = phen_names 
colnames(np_MZm) = phen_names; rownames(np_MZm) = phen_names
colnames(np_DZw) = phen_names; rownames(np_DZw) = phen_names 
colnames(np_DZm) = phen_names; rownames(np_DZm) = phen_names
colnames(np_DZos)= phen_names; rownames(np_DZos) = phen_names


cor_np_list = list(
round(cor_MZw,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(cor_MZm,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(cor_DZw,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(cor_DZm,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(cor_DZos,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(np_MZw,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(np_MZm,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(np_DZw,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(np_DZm,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen"),
round(np_DZos,3) %>% as.data.frame() %>%  rownames_to_column(var = "phen")
)

# Create a workbook
wb = createWorkbook()

# Create a list for sheet names
sn = c("MZw","MZm","DZw","DZm","DZos","npMZw","npMZm","npDZw","npDZm","npDZos")

# Add each correlation matrix as a sheet
for (i in 1:length(cor_np_list)) {
  sheet_name <- sn[i]
  addWorksheet(wb, sheet_name)  # Add a new sheet
  writeData(wb, sheet = sheet_name, x = cor_np_list[[i]]) 
}

# save Supplementary file
saveWorkbook(wb,sprintf("%s/%s/04_supplementary_data_3.xlsx", wdOA,wdOA_sfile), overwrite = T)