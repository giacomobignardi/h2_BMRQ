#SATURATED MODEL####
# structural equation model specification for 5 groups and one covariate:
# 1:monozygotic female
# 2:monozygotic male
# 3:dizygotic female
# 4:dizygotic male
# 5:dizygotic opposite-sex
# each N model is recursively constrained

SAT_model_5v <-"
      #means
    P1_T1~c(mean_f_mz_P1_T1, mean_m_mz_P1_T1, mean_f_dz_P1_T1, mean_m_dz_P1_T1,mean_f_dos_P1_T1)*1 + b_P1*age
    P1_T2~c(mean_f_mz_P1_T2, mean_m_mz_P1_T2, mean_f_dz_P1_T2, mean_m_dz_P1_T2,mean_m_dos_P1_T2)*1 + b_P1*age
    P2_T1~c(mean_f_mz_P2_T1, mean_m_mz_P2_T1, mean_f_dz_P2_T1, mean_m_dz_P2_T1,mean_f_dos_P2_T1)*1 + b_P2*age
    P2_T2~c(mean_f_mz_P2_T2, mean_m_mz_P2_T2, mean_f_dz_P2_T2, mean_m_dz_P2_T2,mean_m_dos_P2_T2)*1 + b_P2*age
    P3_T1~c(mean_f_mz_P3_T1, mean_m_mz_P3_T1, mean_f_dz_P3_T1, mean_m_dz_P3_T1,mean_f_dos_P3_T1)*1 + b_P3*age
    P3_T2~c(mean_f_mz_P3_T2, mean_m_mz_P3_T2, mean_f_dz_P3_T2, mean_m_dz_P3_T2,mean_m_dos_P3_T2)*1 + b_P3*age
    P4_T1~c(mean_f_mz_P4_T1, mean_m_mz_P4_T1, mean_f_dz_P4_T1, mean_m_dz_P4_T1,mean_f_dos_P4_T1)*1 + b_P4*age
    P4_T2~c(mean_f_mz_P4_T2, mean_m_mz_P4_T2, mean_f_dz_P4_T2, mean_m_dz_P4_T2,mean_m_dos_P4_T2)*1 + b_P4*age
    P5_T1~c(mean_f_mz_P5_T1, mean_m_mz_P5_T1, mean_f_dz_P5_T1, mean_m_dz_P5_T1,mean_f_dos_P5_T1)*1 + b_P5*age
    P5_T2~c(mean_f_mz_P5_T2, mean_m_mz_P5_T2, mean_f_dz_P5_T2, mean_m_dz_P5_T2,mean_m_dos_P5_T2)*1 + b_P5*age
    
    #constrain covariate variances to be equal across groups
    age ~~ var_age*age
  
    #variances
    P1_T1~~c(var_f_mz_P1_T1, var_m_mz_P1_T1,var_f_dz_P1_T1, var_m_dz_P1_T1,var_f_dos_P1_T1)*P1_T1
    P1_T2~~c(var_f_mz_P1_T2, var_m_mz_P1_T2,var_f_dz_P1_T2, var_m_dz_P1_T2,var_m_dos_P1_T2)*P1_T2
    P2_T1~~c(var_f_mz_P2_T1, var_m_mz_P2_T1,var_f_dz_P2_T1, var_m_dz_P2_T1,var_f_dos_P2_T1)*P2_T1
    P2_T2~~c(var_f_mz_P2_T2, var_m_mz_P2_T2,var_f_dz_P2_T2, var_m_dz_P2_T2,var_m_dos_P2_T2)*P2_T2
    P3_T1~~c(var_f_mz_P3_T1, var_m_mz_P3_T1,var_f_dz_P3_T1, var_m_dz_P3_T1,var_f_dos_P3_T1)*P3_T1
    P3_T2~~c(var_f_mz_P3_T2, var_m_mz_P3_T2,var_f_dz_P3_T2, var_m_dz_P3_T2,var_m_dos_P3_T2)*P3_T2
    P4_T1~~c(var_f_mz_P4_T1, var_m_mz_P4_T1,var_f_dz_P4_T1, var_m_dz_P4_T1,var_f_dos_P4_T1)*P4_T1
    P4_T2~~c(var_f_mz_P4_T2, var_m_mz_P4_T2,var_f_dz_P4_T2, var_m_dz_P4_T2,var_m_dos_P4_T2)*P4_T2
    P5_T1~~c(var_f_mz_P5_T1, var_m_mz_P5_T1,var_f_dz_P5_T1, var_m_dz_P5_T1,var_f_dos_P5_T1)*P5_T1
    P5_T2~~c(var_f_mz_P5_T2, var_m_mz_P5_T2,var_f_dz_P5_T2, var_m_dz_P5_T2,var_m_dos_P5_T2)*P5_T2
    
    
    #covariances
    P1_T1 ~~ 
    c(cov_f_mz_P12_T11, cov_m_mz_P12_T11, cov_f_dz_P12_T11, cov_m_dz_P12_T11, cov_dos_P12_T11)*P2_T1 + 
    c(cov_f_mz_P13_T11, cov_m_mz_P13_T11, cov_f_dz_P13_T11, cov_m_dz_P13_T11, cov_dos_P13_T11)*P3_T1 + 
    c(cov_f_mz_P14_T11, cov_m_mz_P14_T11, cov_f_dz_P14_T11, cov_m_dz_P14_T11, cov_dos_P14_T11)*P4_T1 +
    c(cov_f_mz_P15_T11, cov_m_mz_P15_T11, cov_f_dz_P15_T11, cov_m_dz_P15_T11, cov_dos_P15_T11)*P5_T1 +
    c(cov_f_mz_P11_T12, cov_m_mz_P11_T12, cov_f_dz_P11_T12, cov_m_dz_P11_T12, cov_dos_P11_T12)*P1_T2 + 
    c(cov_f_mz_P12_T12, cov_m_mz_P12_T12, cov_f_dz_P12_T12, cov_m_dz_P12_T12, cov_dos_P12_T12)*P2_T2 + 
    c(cov_f_mz_P13_T12, cov_m_mz_P13_T12, cov_f_dz_P13_T12, cov_m_dz_P13_T12, cov_dos_P13_T12)*P3_T2 + 
    c(cov_f_mz_P14_T12, cov_m_mz_P14_T12, cov_f_dz_P14_T12, cov_m_dz_P14_T12, cov_dos_P14_T12)*P4_T2 +
    c(cov_f_mz_P15_T12, cov_m_mz_P15_T12, cov_f_dz_P15_T12, cov_m_dz_P15_T12, cov_dos_P15_T12)*P5_T2 

    P2_T1 ~~ 
    c(cov_f_mz_P23_T11, cov_m_mz_P23_T11, cov_f_dz_P23_T11, cov_m_dz_P23_T11, cov_dos_P23_T11)*P3_T1 + 
    c(cov_f_mz_P24_T11, cov_m_mz_P24_T11, cov_f_dz_P24_T11, cov_m_dz_P24_T11, cov_dos_P24_T11)*P4_T1 +
    c(cov_f_mz_P25_T11, cov_m_mz_P25_T11, cov_f_dz_P25_T11, cov_m_dz_P25_T11, cov_dos_P25_T11)*P5_T1 +
    c(cov_f_mz_P21_T12, cov_m_mz_P21_T12, cov_f_dz_P21_T12, cov_m_dz_P21_T12, cov_dos_P21_T12)*P1_T2 + 
    c(cov_f_mz_P22_T12, cov_m_mz_P22_T12, cov_f_dz_P22_T12, cov_m_dz_P22_T12, cov_dos_P22_T12)*P2_T2 + 
    c(cov_f_mz_P23_T12, cov_m_mz_P23_T12, cov_f_dz_P23_T12, cov_m_dz_P23_T12, cov_dos_P23_T12)*P3_T2 + 
    c(cov_f_mz_P24_T12, cov_m_mz_P24_T12, cov_f_dz_P24_T12, cov_m_dz_P24_T12, cov_dos_P24_T12)*P4_T2 +
    c(cov_f_mz_P25_T12, cov_m_mz_P25_T12, cov_f_dz_P25_T12, cov_m_dz_P25_T12, cov_dos_P25_T12)*P5_T2
    
    P3_T1 ~~ 
    c(cov_f_mz_P34_T11, cov_m_mz_P34_T11, cov_f_dz_P34_T11, cov_m_dz_P34_T11, cov_dos_P34_T11)*P4_T1 +
    c(cov_f_mz_P35_T11, cov_m_mz_P35_T11, cov_f_dz_P35_T11, cov_m_dz_P35_T11, cov_dos_P35_T11)*P5_T1 + 
    c(cov_f_mz_P31_T12, cov_m_mz_P31_T12, cov_f_dz_P31_T12, cov_m_dz_P31_T12, cov_dos_P31_T12)*P1_T2 + 
    c(cov_f_mz_P32_T12, cov_m_mz_P32_T12, cov_f_dz_P32_T12, cov_m_dz_P32_T12, cov_dos_P32_T12)*P2_T2 + 
    c(cov_f_mz_P33_T12, cov_m_mz_P33_T12, cov_f_dz_P33_T12, cov_m_dz_P33_T12, cov_dos_P33_T12)*P3_T2 + 
    c(cov_f_mz_P34_T12, cov_m_mz_P34_T12, cov_f_dz_P34_T12, cov_m_dz_P34_T12, cov_dos_P34_T12)*P4_T2 +
    c(cov_f_mz_P35_T12, cov_m_mz_P35_T12, cov_f_dz_P35_T12, cov_m_dz_P35_T12, cov_dos_P35_T12)*P5_T2 
    
    P4_T1 ~~ 
    c(cov_f_mz_P45_T11, cov_m_mz_P45_T11, cov_f_dz_P45_T11, cov_m_dz_P45_T11, cov_dos_P45_T11)*P5_T1 +
    c(cov_f_mz_P41_T12, cov_m_mz_P41_T12, cov_f_dz_P41_T12, cov_m_dz_P41_T12, cov_dos_P41_T12)*P1_T2 + 
    c(cov_f_mz_P42_T12, cov_m_mz_P42_T12, cov_f_dz_P42_T12, cov_m_dz_P42_T12, cov_dos_P42_T12)*P2_T2 + 
    c(cov_f_mz_P43_T12, cov_m_mz_P43_T12, cov_f_dz_P43_T12, cov_m_dz_P43_T12, cov_dos_P43_T12)*P3_T2 + 
    c(cov_f_mz_P44_T12, cov_m_mz_P44_T12, cov_f_dz_P44_T12, cov_m_dz_P44_T12, cov_dos_P44_T12)*P4_T2 +
    c(cov_f_mz_P45_T12, cov_m_mz_P45_T12, cov_f_dz_P45_T12, cov_m_dz_P45_T12, cov_dos_P45_T12)*P5_T2 
    
    P5_T1 ~~ 
    c(cov_f_mz_P51_T12, cov_m_mz_P51_T12, cov_f_dz_P51_T12, cov_m_dz_P51_T12, cov_dos_P51_T12)*P1_T2 + 
    c(cov_f_mz_P52_T12, cov_m_mz_P52_T12, cov_f_dz_P52_T12, cov_m_dz_P52_T12, cov_dos_P52_T12)*P2_T2 + 
    c(cov_f_mz_P53_T12, cov_m_mz_P53_T12, cov_f_dz_P53_T12, cov_m_dz_P53_T12, cov_dos_P53_T12)*P3_T2 + 
    c(cov_f_mz_P54_T12, cov_m_mz_P54_T12, cov_f_dz_P54_T12, cov_m_dz_P54_T12, cov_dos_P54_T12)*P4_T2 +
    c(cov_f_mz_P55_T12, cov_m_mz_P55_T12, cov_f_dz_P55_T12, cov_m_dz_P55_T12, cov_dos_P55_T12)*P5_T2 
    
    P1_T2 ~~ 
    c(cov_f_mz_P12_22, cov_m_mz_P12_22, cov_f_dz_P12_22, cov_m_dz_P12_22, cov_dos_P12_22)*P2_T2 + 
    c(cov_f_mz_P13_22, cov_m_mz_P13_22, cov_f_dz_P13_22, cov_m_dz_P13_22, cov_dos_P13_22)*P3_T2 + 
    c(cov_f_mz_P14_22, cov_m_mz_P14_22, cov_f_dz_P14_22, cov_m_dz_P14_22, cov_dos_P14_22)*P4_T2 +
    c(cov_f_mz_P15_22, cov_m_mz_P15_22, cov_f_dz_P15_22, cov_m_dz_P15_22, cov_dos_P15_22)*P5_T2

    P2_T2 ~~ 
    c(cov_f_mz_P23_22, cov_m_mz_P23_22, cov_f_dz_P23_22, cov_m_dz_P23_22, cov_dos_P23_22)*P3_T2 + 
    c(cov_f_mz_P24_22, cov_m_mz_P24_22, cov_f_dz_P24_22, cov_m_dz_P24_22, cov_dos_P24_22)*P4_T2 +
    c(cov_f_mz_P25_22, cov_m_mz_P25_22, cov_f_dz_P25_22, cov_m_dz_P25_22, cov_dos_P25_22)*P5_T2 
    
    P3_T2 ~~ c(cov_f_mz_P34_22, cov_m_mz_P34_22,cov_f_dz_P34_22, cov_m_dz_P34_22, cov_dos_P34_22)*P4_T2 +
    c(cov_f_mz_P35_22, cov_m_mz_P35_22,cov_f_dz_P35_22, cov_m_dz_P35_22, cov_dos_P35_22)*P5_T2
    
    P4_T2 ~~ c(cov_f_mz_P45_22, cov_m_mz_P45_22,cov_f_dz_P45_22, cov_m_dz_P45_22, cov_dos_P45_22)*P5_T2
    
    #tidy phenotypic correlations
    rp1_f_mz := cov_f_mz_P11_T12 / sqrt(var_f_mz_P1_T1*var_f_mz_P1_T2)
    rp1_m_mz := cov_m_mz_P11_T12 / sqrt(var_m_mz_P1_T1*var_m_mz_P1_T2)
    rp1_f_dz := cov_f_dz_P11_T12 / sqrt(var_f_dz_P1_T1*var_f_dz_P1_T2)
    rp1_m_dz := cov_m_dz_P11_T12 / sqrt(var_m_dz_P1_T1*var_m_dz_P1_T2)
    rp1_dos := cov_dos_P11_T12 / sqrt(var_f_dos_P1_T1*var_m_dos_P1_T2)
    
    rp2_f_mz := cov_f_mz_P22_T12 / sqrt(var_f_mz_P2_T1*var_f_mz_P2_T2)
    rp2_m_mz := cov_m_mz_P22_T12 / sqrt(var_m_mz_P2_T1*var_m_mz_P2_T2)
    rp2_f_dz := cov_f_dz_P22_T12 / sqrt(var_f_dz_P2_T1*var_f_dz_P2_T2)
    rp2_m_dz := cov_m_dz_P22_T12 / sqrt(var_m_dz_P2_T1*var_m_dz_P2_T2)
    rp2_dos := cov_dos_P22_T12 / sqrt(var_f_dos_P2_T1*var_m_dos_P2_T2)
    
    rp3_f_mz := cov_f_mz_P33_T12 / sqrt(var_f_mz_P3_T1*var_f_mz_P3_T2)
    rp3_m_mz := cov_m_mz_P33_T12 / sqrt(var_m_mz_P3_T1*var_m_mz_P3_T2)
    rp3_f_dz := cov_f_dz_P33_T12 / sqrt(var_f_dz_P3_T1*var_f_dz_P3_T2)
    rp3_m_dz := cov_m_dz_P33_T12 / sqrt(var_m_dz_P3_T1*var_m_dz_P3_T2)
    rp3_dos := cov_dos_P33_T12 / sqrt(var_f_dos_P3_T1*var_m_dos_P3_T2)
    
    rp4_f_mz := cov_f_mz_P44_T12 / sqrt(var_f_mz_P4_T1*var_f_mz_P4_T2)
    rp4_m_mz := cov_m_mz_P44_T12 / sqrt(var_m_mz_P4_T1*var_m_mz_P4_T2)
    rp4_f_dz := cov_f_dz_P44_T12 / sqrt(var_f_dz_P4_T1*var_f_dz_P4_T2)
    rp4_m_dz := cov_m_dz_P44_T12 / sqrt(var_m_dz_P4_T1*var_m_dz_P4_T2)
    rp4_dos := cov_dos_P44_T12 / sqrt(var_f_dos_P4_T1*var_m_dos_P4_T2)
    
    rp5_f_mz := cov_f_mz_P55_T12 / sqrt(var_f_mz_P5_T1*var_f_mz_P5_T2)
    rp5_m_mz := cov_m_mz_P55_T12 / sqrt(var_m_mz_P5_T1*var_m_mz_P5_T2)
    rp5_f_dz := cov_f_dz_P55_T12 / sqrt(var_f_dz_P5_T1*var_f_dz_P5_T2)
    rp5_m_dz := cov_m_dz_P55_T12 / sqrt(var_m_dz_P5_T1*var_m_dz_P5_T2)
    rp5_dos := cov_dos_P55_T12 / sqrt(var_f_dos_P5_T1*var_m_dos_P5_T2)
 "

# use this model to extract phenotypic correlations and cross-trait cross-twin correlations
SAT_model_5v_const <-"
      #means
    P1_T1~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1,mean_f_P1)*1 + b_P1*age
    P1_T2~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1,mean_m_P1)*1 + b_P1*age
    P2_T1~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2,mean_f_P2)*1 + b_P2*age
    P2_T2~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2,mean_m_P2)*1 + b_P2*age
    P3_T1~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3,mean_f_P3)*1 + b_P3*age
    P3_T2~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3,mean_m_P3)*1 + b_P3*age
    P4_T1~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4,mean_f_P4)*1 + b_P4*age
    P4_T2~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4,mean_m_P4)*1 + b_P4*age
    P5_T1~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5,mean_f_P5)*1 + b_P5*age
    P5_T2~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5,mean_m_P5)*1 + b_P5*age
    
    #constrain covariate variances to be equal across groups
    age ~~ var_age*age
  
    #variances
    P1_T1~~c(var_P1)*P1_T1
    P1_T2~~c(var_P1)*P1_T2
    P2_T1~~c(var_P2)*P2_T1
    P2_T2~~c(var_P2)*P2_T2
    P3_T1~~c(var_P3)*P3_T1
    P3_T2~~c(var_P3)*P3_T2
    P4_T1~~c(var_P4)*P4_T1
    P4_T2~~c(var_P4)*P4_T2
    P5_T1~~c(var_P5)*P5_T1
    P5_T2~~c(var_P5)*P5_T2
    
    
    #covariances
    P1_T1 ~~ 
    c(cov_mz_P12, cov_mz_P12, cov_dz_P12, cov_dz_P12, cov_dz_P12)*P2_T1 + 
    c(cov_mz_P13, cov_mz_P13, cov_dz_P13, cov_dz_P13, cov_dz_P13)*P3_T1 + 
    c(cov_mz_P14, cov_mz_P14, cov_dz_P14, cov_dz_P14, cov_dz_P14)*P4_T1 +
    c(cov_mz_P15, cov_mz_P15, cov_dz_P15, cov_dz_P15, cov_dz_P15)*P5_T1 +
    c(cov_mz_P11_T12, cov_mz_P11_T12, cov_dz_P11_T12, cov_dz_P11_T12, cov_dz_P11_T12)*P1_T2 + 
    c(cov_mz_P12_T12, cov_mz_P12_T12, cov_dz_P12_T12, cov_dz_P12_T12, cov_dz_P12_T12)*P2_T2 + 
    c(cov_mz_P13_T12, cov_mz_P13_T12, cov_dz_P13_T12, cov_dz_P13_T12, cov_dz_P13_T12)*P3_T2 + 
    c(cov_mz_P14_T12, cov_mz_P14_T12, cov_dz_P14_T12, cov_dz_P14_T12, cov_dz_P14_T12)*P4_T2 +
    c(cov_mz_P15_T12, cov_mz_P15_T12, cov_dz_P15_T12, cov_dz_P15_T12, cov_dz_P15_T12)*P5_T2 

    P2_T1 ~~ 
    c(cov_mz_P23, cov_mz_P23, cov_dz_P23, cov_dz_P23, cov_dz_P23)*P3_T1 + 
    c(cov_mz_P24, cov_mz_P24, cov_dz_P24, cov_dz_P24, cov_dz_P24)*P4_T1 +
    c(cov_mz_P25, cov_mz_P25, cov_dz_P25, cov_dz_P25, cov_dz_P25)*P5_T1 +
    c(cov_mz_P12_T12, cov_mz_P12_T12, cov_dz_P12_T12, cov_dz_P12_T12, cov_dz_P12_T12)*P1_T2 + 
    c(cov_mz_P22_T12, cov_mz_P22_T12, cov_dz_P22_T12, cov_dz_P22_T12, cov_dz_P22_T12)*P2_T2 + 
    c(cov_mz_P23_T12, cov_mz_P23_T12, cov_dz_P23_T12, cov_dz_P23_T12, cov_dz_P23_T12)*P3_T2 + 
    c(cov_mz_P24_T12, cov_mz_P24_T12, cov_dz_P24_T12, cov_dz_P24_T12, cov_dz_P24_T12)*P4_T2 +
    c(cov_mz_P25_T12, cov_mz_P25_T12, cov_dz_P25_T12, cov_dz_P25_T12, cov_dz_P25_T12)*P5_T2
    
    P3_T1 ~~ 
    c(cov_mz_P34, cov_mz_P34, cov_dz_P34, cov_dz_P34, cov_dz_P34)*P4_T1 +
    c(cov_mz_P35, cov_mz_P35, cov_dz_P35, cov_dz_P35, cov_dz_P35)*P5_T1 + 
    c(cov_mz_P13_T12, cov_mz_P13_T12, cov_dz_P13_T12, cov_dz_P13_T12, cov_dz_P13_T12)*P1_T2 + 
    c(cov_mz_P23_T12, cov_mz_P23_T12, cov_dz_P23_T12, cov_dz_P23_T12, cov_dz_P23_T12)*P2_T2 + 
    c(cov_mz_P33_T12, cov_mz_P33_T12, cov_dz_P33_T12, cov_dz_P33_T12, cov_dz_P33_T12)*P3_T2 + 
    c(cov_mz_P34_T12, cov_mz_P34_T12, cov_dz_P34_T12, cov_dz_P34_T12, cov_dz_P34_T12)*P4_T2 +
    c(cov_mz_P35_T12, cov_mz_P35_T12, cov_dz_P35_T12, cov_dz_P35_T12, cov_dz_P35_T12)*P5_T2 
    
    P4_T1 ~~ 
    c(cov_mz_P45, cov_mz_P45, cov_dz_P45, cov_dz_P45, cov_dz_P45)*P5_T1 +
    c(cov_mz_P14_T12, cov_mz_P14_T12, cov_dz_P14_T12, cov_dz_P14_T12, cov_dz_P14_T12)*P1_T2 + 
    c(cov_mz_P24_T12, cov_mz_P24_T12, cov_dz_P24_T12, cov_dz_P24_T12, cov_dz_P24_T12)*P2_T2 + 
    c(cov_mz_P34_T12, cov_mz_P34_T12, cov_dz_P34_T12, cov_dz_P34_T12, cov_dz_P34_T12)*P3_T2 + 
    c(cov_mz_P44_T12, cov_mz_P44_T12, cov_dz_P44_T12, cov_dz_P44_T12, cov_dz_P44_T12)*P4_T2 +
    c(cov_mz_P45_T12, cov_mz_P45_T12, cov_dz_P45_T12, cov_dz_P45_T12, cov_dz_P45_T12)*P5_T2 
    
    P5_T1 ~~ 
    c(cov_mz_P15_T12, cov_mz_P15_T12, cov_dz_P15_T12, cov_dz_P15_T12, cov_dz_P15_T12)*P1_T2 + 
    c(cov_mz_P25_T12, cov_mz_P25_T12, cov_dz_P25_T12, cov_dz_P25_T12, cov_dz_P25_T12)*P2_T2 + 
    c(cov_mz_P35_T12, cov_mz_P35_T12, cov_dz_P35_T12, cov_dz_P35_T12, cov_dz_P35_T12)*P3_T2 + 
    c(cov_mz_P45_T12, cov_mz_P45_T12, cov_dz_P45_T12, cov_dz_P45_T12, cov_dz_P45_T12)*P4_T2 +
    c(cov_mz_P55_T12, cov_mz_P55_T12, cov_dz_P55_T12, cov_dz_P55_T12, cov_dz_P55_T12)*P5_T2 
    
    P1_T2 ~~ 
    c(cov_mz_P12, cov_mz_P12, cov_dz_P12, cov_dz_P12, cov_dz_P12)*P2_T2 + 
    c(cov_mz_P13, cov_mz_P13, cov_dz_P13, cov_dz_P13, cov_dz_P13)*P3_T2 + 
    c(cov_mz_P14, cov_mz_P14, cov_dz_P14, cov_dz_P14, cov_dz_P14)*P4_T2 +
    c(cov_mz_P15, cov_mz_P15, cov_dz_P15, cov_dz_P15, cov_dz_P15)*P5_T2

    P2_T2 ~~ 
    c(cov_mz_P23, cov_mz_P23, cov_dz_P23, cov_dz_P23, cov_dz_P23)*P3_T2 + 
    c(cov_mz_P24, cov_mz_P24, cov_dz_P24, cov_dz_P24, cov_dz_P24)*P4_T2 +
    c(cov_mz_P25, cov_mz_P25, cov_dz_P25, cov_dz_P25, cov_dz_P25)*P5_T2 
    
    P3_T2 ~~ c(cov_mz_P34, cov_mz_P34,cov_dz_P34, cov_dz_P34, cov_dz_P34)*P4_T2 +
    c(cov_mz_P35, cov_mz_P35,cov_dz_P35, cov_dz_P35, cov_dz_P35)*P5_T2
    
    P4_T2 ~~ c(cov_mz_P45, cov_mz_P45,cov_dz_P45, cov_dz_P45, cov_dz_P45)*P5_T2
    
    #tidy twin phenotypic correlations
    rp1_mz := cov_mz_P11_T12 / sqrt(var_P1*var_P1)
    rp1_dz := cov_dz_P11_T12 / sqrt(var_P1*var_P1)
    
    rp2_mz := cov_mz_P22_T12 / sqrt(var_P2*var_P2)
    rp2_dz := cov_dz_P22_T12 / sqrt(var_P2*var_P2)
    
    rp3_mz := cov_mz_P33_T12 / sqrt(var_P3*var_P3)
    rp3_dz := cov_dz_P33_T12 / sqrt(var_P3*var_P3)
    
    rp4_mz := cov_mz_P44_T12 / sqrt(var_P4*var_P4)
    rp4_dz := cov_dz_P44_T12 / sqrt(var_P4*var_P4)
    
    rp5_mz := cov_mz_P55_T12 / sqrt(var_P5*var_P5)
    rp5_dz := cov_dz_P55_T12 / sqrt(var_P5*var_P5)
    
 "