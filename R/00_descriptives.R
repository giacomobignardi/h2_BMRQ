#Author: Giacomo Bignardi
#Adapted from: NA
#
#
#
#Description: Prepare data for later analysis
#TIDY: Tidy initial dataframe
#PHENOTYPES: prepare phenotype (bmrq & facets, reward responsiveness, music discrimination)
#Program: 00 ------------------------------------------------------------------------------------------------------------------------------
.libPaths('W:\\XXX\04_packages')

# load packages
library(tidyverse)
library(farver)
library(patchwork)

# clean working environment 
rm(list = ls())

# set Open Access working directories
wdOA = getwd()
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs"

# set not Open Access working directories
wdNOA_Data = "01_data"
wdOA_ImageOutput = "05_images"

# load dataFrames:
twin_BMRQ  = read_csv(sprintf("%s/%s/XXX.csv", wdOA,wdNOA_Data))

# sanity check that there are no duplicated families
unique(table(twin_BMRQ$pairnr))
nrow(twin_BMRQ)

# set seed for reproducibility
set.seed(123)

# TIDY####
# convert null to NA
twin_BMRQ = twin_BMRQ %>% 
  mutate(rhythm_score = as.numeric(rhythm_score),
         melody_d_score = as.numeric(melody_d_score),
         pitch_score = as.numeric(pitch_score)
  )

# remove twins with uncertain zygosity
twin_BMRQ = twin_BMRQ %>% filter(bestzyg != 3)

# check that there are no missing SEX
sum(is.na(twin_BMRQ$sex))

# check that there are no missing zygosity
sum(is.na(twin_BMRQ$bestzyg))
table(twin_BMRQ$bestzyg)

# check that there are no mislabeled zygosity
twin_BMRQ = twin_BMRQ %>% 
  mutate(unique_pairzyg = paste0(as.character(pairnr), as.character(bestzyg))) %>%
  mutate(pair = duplicated(unique_pairzyg)) %>%
  mutate(unique_pair = duplicated(pairnr))
sum(!(twin_BMRQ$unique_pair == twin_BMRQ$unique_pair))

# recode twin zygosities (1 MZfemale; 2 MZmale; 3 DZfemale, 4 DZmale; 5 DZoppositesex)
twin_BMRQ = twin_BMRQ %>% 
  mutate(zyg = ifelse(bestzyg == 1 & sex == 1,2, #MZmale 
                      ifelse(bestzyg == 1 & sex == 2,1, #MZfemale
                             ifelse(bestzyg == 2 & sex == 1,4, #DZmale
                                    ifelse(bestzyg == 2 & sex == 2,3, #DZfemale
                                           ifelse(bestzyg == 4, 5, NA #DOS 
                                             )))))) %>% 
  rename(fam = pairnr)
table(twin_BMRQ$zyg)

# recode age
twin_BMRQ$age = round(as.numeric(
  difftime(
    paste0(substr(twin_BMRQ$SURVEY_DATE,1,7),"-01"), #get month and year of completion of survey
    paste0(substr(twin_BMRQ$BIRTH_YYYYMM, 1, 4) ,"-",substr(twin_BMRQ$BIRTH_YYYYMM, 5, 6),"-01"), #get month and year of age
    unit="weeks") / 52.25 #divide by average number of weeks per year
),0)


# PHENOTYPES####
# Wave 1 measures
# Rhythm, melody, and pitch discrimination score; Ullen et al., (2014) https://www.sciencedirect.com/science/article/pii/S0191886914000841
# rhythm_score
# melody_d_score
# pitch_score

# Wave 2 measures
# BMRQ answer options (Mas-Herrero et al,. 2013): english; sweden
# [1] - Completely disagree; Håller inte alls med
# [2] - disagree; Håller delvis inte med
# [3] - Neither agree nor disagree; Varken eller
# [4] - agree; Håller delvis med
# [5] Completely agree; Håller helt med
# BIS-BAS answer options (Carver & White, 1994); https://local.psy.miami.edu/people/faculty/ccarver/availbale-self-report-instruments/bisbas-scales/ ): english; sweden
# [1] - very true for me; Mycket sant för mig
# [2] - somewhat true for me; Något sant för mig
# [3] - somewhat false for me; Något falskt för mig
# [4] - very false for me; Mycket falskt för mig
# NOTE THAT BFI DATA ARE USED FOR ANOTHER PROJECT!
# BFI answer options (Soto & John, 2017). https://www.sciencedirect.com/science/article/pii/S0092656616301325 (note that the scoring is different in the STAGE version) english; sweden
# [0] - Completely disagree; Håller inte alls med
# [1] - disagree; Håller delvis inte med
# [2] - Neither agree nor disagree; Varken eller
# [3] - agree; Håller delvis med
# [4] - Completely agree; Håller helt med

# reverse score for reverse items
twin_BMRQ = twin_BMRQ %>% mutate(
                                 #BMRQ reversed items
                                 #Musical Seeking
                                 BMRS_2_r = 6-BMRS_2,
                                 #Sensory Motor
                                 BMRS_5_r = 6-BMRS_5,
                          
                                 #note that for the BIS/BAS the order of answer is already reversed 
                                 #BAS Reward Responsiveness
                                 BISBAS_4_r = 5-BISBAS_4,
                                 BISBAS_7_r = 5-BISBAS_7,
                                 BISBAS_14_r = 5-BISBAS_14,
                                 BISBAS_18_r = 5-BISBAS_18,
                                 BISBAS_23_r = 5-BISBAS_23,
                                 )

# create sum score and parcels
twin_BMRQ = twin_BMRQ %>% mutate(#overall score
  BMRS_total = 
    BMRS_1 +
    BMRS_2_r +
    BMRS_3 +
    BMRS_4 +
    BMRS_5_r +
    BMRS_6 +
    BMRS_7 +
    BMRS_8 +
    BMRS_9 +
    BMRS_10 +
    BMRS_11 +
    BMRS_12 +
    BMRS_13 +
    BMRS_14 +
    BMRS_15 +
    BMRS_16 +
    BMRS_17 +
    BMRS_18 +
    BMRS_19 +
    BMRS_20,
  #Musical Seeking score
  BMRS_MS = 
    BMRS_11 +
    BMRS_7 +
    BMRS_17 +
    BMRS_2_r,
  #Emotion Evocation score
  BMRS_EE = 
    BMRS_18 +
    BMRS_12 +
    BMRS_8 +
    BMRS_3,
  #Mood Regulation score
  BMRS_MR = 
    BMRS_14 +
    BMRS_9 +
    BMRS_19 +
    BMRS_4,
  #Senosry Motor score
  BMRS_SM = 
    BMRS_10 +
    BMRS_20+
    BMRS_15 +
    BMRS_5_r,
  #Social Reward score
  BMRS_SR = 
    BMRS_13 +
    BMRS_1 +
    BMRS_6 +
    BMRS_16,
  #BAS Reward Responsiveness
  BAS_RR = 
    BISBAS_4_r + 
    BISBAS_7_r + 
    BISBAS_14_r + 
    BISBAS_18_r + 
    BISBAS_23_r
)

# remove individuals for which the BMRQ is missing
twin_BMRQ = twin_BMRQ %>% filter(!(is.na(BMRS_total)))

# randomize twin order by randomizing the pair number 
paired_index = twin_BMRQ %>% filter(duplicated(fam)) %>% pull(fam)
# sanity check 1
sum(table(twin_BMRQ$tvab))
twin_BMRQ = twin_BMRQ%>% 
  mutate(paired = ifelse(fam %in% paired_index, 1,0)) %>% 
  group_by(fam) %>% 
  #randomly assign order to twin pair with a non missing pair
  mutate(tvab = ifelse(paired == 0, tvab, sample(tvab))) %>% 
  #order n pair for DOS: if male then pair 2, if female than pair 1
  mutate(tvab = ifelse(zyg ==5 & sex == 2 & tvab == 2, 1,
                       ifelse(zyg ==5 & sex == 1  & tvab == 1, 2,tvab))) %>%
  ungroup()
# sanity check 2
sum(table(twin_BMRQ$tvab))

# WIDE FORMAT####
# prepare wide format
twin_BMRQ_wide = twin_BMRQ %>% select(
  fam, 
  tvab,#to later wide the dataframe
  zyg,
  age,
  starts_with("BMRS"),
  BAS_RR,
  music_discrimination,
  rhythm_score,
  melody_d_score,
  pitch_score) %>% 
  pivot_wider(names_from = tvab, values_from = age:pitch_score)

# inspect age (need to be exactly the same)
twin_BMRQ_wide %>% ggplot(aes(age_1,age_2)) + geom_point()

# remove age2 as it is redundant
twin_BMRQ_wide = twin_BMRQ_wide %>% 
  mutate(
    age_1 =  ifelse(is.na(age_1) , age_2, age_1),
    age_2 = ifelse(is.na(age_2), age_1, age_2)
  ) %>% 
  mutate (age = (age_1 + age_2) /2) %>% 
  select(-c(age_1, age_2))

# DESCRIPTIVES#### 
# AGE
# age of the sample
twin_BMRQ %>% rstatix::get_summary_stats(age)%>% select(min,max,mean,sd)
# mean of the BMRQ 
twin_BMRQ %>% rstatix::get_summary_stats(BMRS_total)
# skewness
psych::skew(twin_BMRQ$BMRS_total)

# AGE
# age of the sample
twin_BMRQ %>% rstatix::get_summary_stats(age)%>% select(min,max,mean,sd)
twin_BMRQ %>% rstatix::get_summary_stats(BMRS_total)

# CFA#### 
# CFA to confirm that BMRQ sum score is appropriate in swedish sample
MR_T1_mod =
" MRS =~  NA*BMRS_EE_1 + BMRS_MS_1 + BMRS_MR_1 + BMRS_SM_1 + BMRS_SR_1
  MRS ~~ 1*MRS
"
MR_T1_fit = lavaan::cfa(MR_T1_mod, data = twin_BMRQ_wide)
lavaan::fitmeasures(MR_T1_fit, c("cfi", "srmr"))

MR_T2_mod =
  " MRS =~  NA*BMRS_EE_2 + BMRS_MS_2 + BMRS_MR_2 + BMRS_SM_2 + BMRS_SR_2
  MRS ~~ 1*MRS
"
MR_T2_fit = lavaan::cfa(MR_T2_mod, data = twin_BMRQ_wide)
lavaan::fitmeasures(MR_T2_fit, c("cfi", "srmr"))

# N individuals with full BMRQ data
twin_BMRQ %>% filter(!is.na(BMRS_total)) %>%  reframe(table(zyg))
sum(twin_BMRQ %>% filter(!is.na(BMRS_total)) %>%  reframe(table(zyg)))

# N individuals with full BAS_RR data
twin_BMRQ %>% filter(!is.na(BAS_RR)) %>%  reframe(table(zyg))
sum(twin_BMRQ %>% filter(!is.na(BAS_RR)) %>%  reframe(table(zyg)))

# N individuals with full discrimination scores
twin_BMRQ %>% filter(!is.na(music_discrimination)) %>%  reframe(table(zyg))
sum(twin_BMRQ %>% filter(!is.na(music_discrimination)) %>%  reframe(table(zyg)))

# frequency of MZ and DZ pairs with full BMRQ data
twin_BMRQ_nPairs_1 = twin_BMRQ_wide %>% reframe(table(zyg))
sum(twin_BMRQ_nPairs_1$`table(zyg)`)

# number of complete pairs
twin_BMRQ_wide %>% filter(!(is.na(BMRS_total_1) | is.na(BMRS_total_2))) %>%  reframe(table(zyg))
sum(twin_BMRQ_wide %>% filter(!(is.na(BMRS_total_1) | is.na(BMRS_total_2))) %>%  reframe(table(zyg)))

# number of  pairs with BAS_RR_1
twin_BMRQ_wide %>% filter(!(is.na(BAS_RR_1) & is.na(BAS_RR_2))) %>%  reframe(table(zyg))
sum(twin_BMRQ_wide %>% filter(!(is.na(BAS_RR_1) & is.na(BAS_RR_2))) %>%  reframe(table(zyg)))

# number of complete pairs with BAS_RR_1
twin_BMRQ_wide %>% 
  filter(!(is.na(BMRS_total_1) | is.na(BMRS_total_2))) %>% 
  filter(!(is.na(BAS_RR_1) | is.na(BAS_RR_2))) %>%
  reframe(table(zyg))
sum(twin_BMRQ_wide %>% 
      filter(!(is.na(BMRS_total_1) | is.na(BMRS_total_2))) %>% 
      filter(!(is.na(BAS_RR_1) | is.na(BAS_RR_2))) %>%
      reframe(table(zyg)))

# number of  pairs with discrimination scores
twin_BMRQ_wide %>% filter(!(is.na(music_discrimination_1) & is.na(music_discrimination_2))) %>%  reframe(table(zyg))
sum(twin_BMRQ_wide %>% filter(!(is.na(music_discrimination_1) & is.na(music_discrimination_2))) %>%  reframe(table(zyg)))

# number of complete pairs with BAS_RR_1 abd discrimination score
twin_BMRQ_wide %>% 
  filter(!(is.na(BMRS_total_1) | is.na(BMRS_total_2))) %>% 
  # filter(!(is.na(BAS_RR_1) | is.na(BAS_RR_2))) %>%
  filter(!(is.na(music_discrimination_1) | is.na(music_discrimination_1))) %>%
  reframe(table(zyg))
sum(twin_BMRQ_wide %>% 
      filter(!(is.na(BMRS_total_1) | is.na(BMRS_total_2))) %>% 
      # filter(!(is.na(BAS_RR_1) | is.na(BAS_RR_2))) %>%
      filter(!(is.na(music_discrimination_1) | is.na(music_discrimination_2))) %>%
      reframe(table(zyg)))



# DISTRIBUTIONS####
# overall
dist_p1_BMRQ = twin_BMRQ_wide %>% 
  ggplot(aes(BMRS_total_1)) + 
  geom_histogram(color = "black", fill = "lightGray", binwidth = 5)+ 
  theme_classic() + 
  ylim(c(0,750)) +
  labs(x = "music reward sensitivity T1")

dist_p2_BMRQ = twin_BMRQ_wide %>% 
  ggplot(aes(BMRS_total_2)) + 
  geom_histogram(color = "black", fill = "lightGray", binwidth = 5)+ 
  theme_classic() + 
  ylim(c(0,750)) +
  labs(title = "") +
  labs(x = "music reward sensitivity T2")

twin_BMRQ_wide %>% summarise(range(BMRS_total_1, na.rm = T))

# SAVE####
# save the wide and cleaned version
write_csv(twin_BMRQ_wide,sprintf("%s/%s/00_twin_BMRQ_wide.csv", wdOA,wdOA_output))

# save the histograms (total sum score)
pdf(sprintf("%s/%s/00_histogram_BMRQ.pdf",wdOA,wdOA_ImageOutput),
    width = 6,
    height =4)
dist_p1_BMRQ|dist_p2_BMRQ
dev.off()