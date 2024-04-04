#SATURATED MODEL####
# structural equation model specification for 5 groups and one covariate:
# 1:monozygotic female
# 2:monozygotic male
# 3:dizygotic female
# 4:dizygotic male
# 5:dizygotic opposite-sex
# each N model is recursively constrained

#full saturated model
SAT_mod =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f_mz_1,mean_m_mz_1,mean_f_dz_1,mean_m_dz_1,mean_f_dos_1)*1 + c(b_f_mz_1,b_m_mz_1,b_f_dz_1,b_m_dz_1,b_f_dos_1)*age
  P_T2 ~ c(mean_f_mz_2,mean_m_mz_2,mean_f_dz_2,mean_m_dz_2,mean_m_dos_2)*1 + c(b_f_mz_2,b_m_mz_2,b_f_dz_2,b_m_dz_2,b_m_dos_2)*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f_mz_1,var_m_mz_1,var_f_dz_1,var_m_dz_1,var_f_dos_1)*P_T1
  P_T2 ~~ c(var_f_mz_2,var_m_mz_2,var_f_dz_2,var_m_dz_2,var_m_dos_2)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
 "


#full saturated model with manifest covariate restriciton
qSAT_mod =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f_mz_1,mean_m_mz_1,mean_f_dz_1,mean_m_dz_1,mean_f_dos_1)*1 + c(b_f_mz_1,b_m_mz_1,b_f_dz_1,b_m_dz_1,b_f_dos_1)*age
  P_T2 ~ c(mean_f_mz_2,mean_m_mz_2,mean_f_dz_2,mean_m_dz_2,mean_m_dos_2)*1 + c(b_f_mz_2,b_m_mz_2,b_f_dz_2,b_m_dz_2,b_m_dos_2)*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f_mz_1,var_m_mz_1,var_f_dz_1,var_m_dz_1,var_f_dos_1)*P_T1
  P_T2 ~~ c(var_f_mz_2,var_m_mz_2,var_f_dz_2,var_m_dz_2,var_m_dos_2)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age  "

#constrained model with same covariate effects and variances 
SAT_mod_cov_0 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f_mz_1,mean_m_mz_1,mean_f_dz_1,mean_m_dz_1,mean_f_dos_1)*1 + b1*age
  P_T2 ~ c(mean_f_mz_2,mean_m_mz_2,mean_f_dz_2,mean_m_dz_2,mean_m_dos_2)*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f_mz_1,var_m_mz_1,var_f_dz_1,var_m_dz_1,var_f_dos_1)*P_T1
  P_T2 ~~ c(var_f_mz_2,var_m_mz_2,var_f_dz_2,var_m_dz_2,var_m_dos_2)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"

#constrained model without covariate
SAT_mod_cov_1 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f_mz_1,mean_m_mz_1,mean_f_dz_1,mean_m_dz_1,mean_f_dos_1)*1  + 0*age
  P_T2 ~ c(mean_f_mz_2,mean_m_mz_2,mean_f_dz_2,mean_m_dz_2,mean_m_dos_2)*1  + 0*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f_mz_1,var_m_mz_1,var_f_dz_1,var_m_dz_1,var_f_dos_1)*P_T1
  P_T2 ~~ c(var_f_mz_2,var_m_mz_2,var_f_dz_2,var_m_dz_2,var_m_dos_2)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"

#constrained model: birth order, implies no birth order effects on means and variances
SAT_mod_bor_2 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f_mz,mean_m_mz,mean_f_dz,mean_m_dz,mean_f_dos)*1 + b1*age
  P_T2 ~ c(mean_f_mz,mean_m_mz,mean_f_dz,mean_m_dz,mean_m_dos)*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f_mz,var_m_mz,var_f_dz,var_m_dz,var_f_dos)*P_T1
  P_T2 ~~ c(var_f_mz,var_m_mz,var_f_dz,var_m_dz,var_m_dos)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"

#constrained model: zygosity, implies no zygosity effects on means and variances
SAT_mod_zyg_3 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
  P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f,var_m,var_f,var_m,var_f)*P_T1
  P_T2 ~~ c(var_f,var_m,var_f,var_m,var_m)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"

#constrained model: sex mean, implies no sex effects on means
SAT_mod_sxm_4 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ mean*1 + b1*age
  P_T2 ~ mean*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ c(var_f,var_m,var_f,var_m,var_f)*P_T1
  P_T2 ~~ c(var_f,var_m,var_f,var_m,var_m)*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"

#constrained model: sex variance, implies no sex effects on variances
SAT_mod_sxv_5 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
  P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ var*P_T1
  P_T2 ~~ var*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_f_mz,cov_m_mz,cov_f_dz,cov_m_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"

#constrained model: covariance sex mz, implies no difference between covariances across same/sex
SAT_mod_qnt_6 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
  P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ var*P_T1
  P_T2 ~~ var*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_mz,cov_mz,cov_dz,cov_dz,cov_dos)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"
#constrained model: covariance sex dz, implies no difference between covariances across opposite-sex and same-sex
SAT_mod_qal_7 =
  "
  #estimate intercepts (means) for MZ and DZ twin 1 and 2 separately
  P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
  P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
  #estimate variances for MZ and DZ twins
  P_T1 ~~ var*P_T1
  P_T2 ~~ var*P_T2
  #estimate covariances within MZ and DZ twin pairs
  P_T1 ~~ c(cov_mz,cov_mz,cov_dz,cov_dz,cov_dz)*P_T2
  #constrain covariate variances to be equal across groups
  age ~~ var_age*age
"
