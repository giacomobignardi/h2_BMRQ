#MULTIVARATE CTD####
# structural equation model specification for 5 groups, 5 phenotypes, and one covariate per phenotype:
# 1:monozygotic female
# 2:monozygotic male
# 3:dizygotic female
# 4:dizygotic male
# 5:dizygotic opposite-sex
#correlated and hybrid independent pathway model, correlated factor solution

cfs_AE_mod <-"
    #means
    P1_T1~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_f_P1)*1 + b_P1*age
    P1_T2~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_m_P1)*1 + b_P1*age
    P2_T1~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_f_P2)*1 + b_P2*age
    P2_T2~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_m_P2)*1 + b_P2*age
    P3_T1~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_f_P3)*1 + b_P3*age
    P3_T2~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_m_P3)*1 + b_P3*age
    P4_T1~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4, mean_f_P4)*1 + b_P4*age
    P4_T2~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4, mean_m_P4)*1 + b_P4*age
    P5_T1~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5, mean_f_P5)*1 + b_P5*age
    P5_T2~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5, mean_m_P5)*1 + b_P5*age
    
    #constrain covariate variances to be equal across groups
    age ~~ var_age*age

    #AE components
    A_P1_T1=~ P1_T1
    A_P1_T2=~ P1_T2
    A_P2_T1=~ P2_T1
    A_P2_T2=~ P2_T2
    A_P3_T1=~ P3_T1
    A_P3_T2=~ P3_T2   
    A_P4_T1=~ P4_T1
    A_P4_T2=~ P4_T2
    A_P5_T1=~ P5_T1
    A_P5_T2=~ P5_T2

    E_P1_T1=~ P1_T1
    E_P1_T2=~ P1_T2
    E_P2_T1=~ P2_T1
    E_P2_T2=~ P2_T2
    E_P3_T1=~ P3_T1
    E_P3_T2=~ P3_T2   
    E_P4_T1=~ P4_T1
    E_P4_T2=~ P4_T2
    E_P5_T1=~ P5_T1
    E_P5_T2=~ P5_T2

    #variances
    A_P1_T1 ~~ varA_P1*A_P1_T1
    A_P1_T2 ~~ varA_P1*A_P1_T2
    A_P2_T1 ~~ varA_P2*A_P2_T1
    A_P2_T2 ~~ varA_P2*A_P2_T2
    A_P3_T1 ~~ varA_P3*A_P3_T1
    A_P3_T2 ~~ varA_P3*A_P3_T2
    A_P4_T1 ~~ varA_P4*A_P4_T1
    A_P4_T2 ~~ varA_P4*A_P4_T2
    A_P5_T1 ~~ varA_P5*A_P5_T1
    A_P5_T2 ~~ varA_P5*A_P5_T2

    E_P1_T1 ~~ varE_P1*E_P1_T1
    E_P1_T2 ~~ varE_P1*E_P1_T2
    E_P2_T1 ~~ varE_P2*E_P2_T1
    E_P2_T2 ~~ varE_P2*E_P2_T2
    E_P3_T1 ~~ varE_P3*E_P3_T1
    E_P3_T2 ~~ varE_P3*E_P3_T2
    E_P4_T1 ~~ varE_P4*E_P4_T1
    E_P4_T2 ~~ varE_P4*E_P4_T2
    E_P5_T1 ~~ varE_P5*E_P5_T1
    E_P5_T2 ~~ varE_P5*E_P5_T2

    #fix remaining variance to 0
    P1_T1~~0*P1_T1
    P1_T2~~0*P1_T2
    P2_T1~~0*P2_T1
    P2_T2~~0*P2_T2
    P3_T1~~0*P3_T1
    P3_T2~~0*P3_T2
    P4_T1~~0*P4_T1
    P4_T2~~0*P4_T2
    P5_T1~~0*P5_T1
    P5_T2~~0*P5_T2
    
    
    #constraints (twin pair covariances, p = phenotype)
    cov_MZ_P1 == varA_P1
    cov_MZ_P2 == varA_P2
    cov_MZ_P3 == varA_P3
    cov_MZ_P4 == varA_P4
    cov_MZ_P5 == varA_P5

    cov_DZ_P1 == .5*varA_P1
    cov_DZ_P2 == .5*varA_P2
    cov_DZ_P3 == .5*varA_P3
    cov_DZ_P4 == .5*varA_P4
    cov_DZ_P5 == .5*varA_P5
    
    #covariances cross-twin wihtin-trait
    A_P1_T1 ~~ c((cov_MZ_P1),(cov_MZ_P1),(cov_DZ_P1),(cov_DZ_P1),(cov_DZ_P1))*A_P1_T2
    A_P2_T1 ~~ c((cov_MZ_P2),(cov_MZ_P2),(cov_DZ_P2),(cov_DZ_P2),(cov_DZ_P2))*A_P2_T2
    A_P3_T1 ~~ c((cov_MZ_P3),(cov_MZ_P3),(cov_DZ_P3),(cov_DZ_P3),(cov_DZ_P3))*A_P3_T2 
    A_P4_T1 ~~ c((cov_MZ_P4),(cov_MZ_P4),(cov_DZ_P4),(cov_DZ_P4),(cov_DZ_P4))*A_P4_T2
    A_P5_T1 ~~ c((cov_MZ_P5),(cov_MZ_P5),(cov_DZ_P5),(cov_DZ_P5),(cov_DZ_P5))*A_P5_T2
    
    E_P1_T1 ~~ 0*E_P1_T2
    E_P2_T1 ~~ 0*E_P2_T2
    E_P3_T1 ~~ 0*E_P3_T2
    E_P4_T1 ~~ 0*E_P4_T2
    E_P5_T1 ~~ 0*E_P5_T2

    #covariances within-twin cross-trait
    A_P1_T1 ~~ (covA_P12)*A_P2_T1 + (covA_P13)*A_P3_T1 + (covA_P14)*A_P4_T1 + (covA_P15)*A_P5_T1
    A_P2_T1 ~~ (covA_P23)*A_P3_T1 + (covA_P24)*A_P4_T1 + (covA_P25)*A_P5_T1
    A_P3_T1 ~~ (covA_P34)*A_P4_T1 + (covA_P35)*A_P5_T1
    A_P4_T1 ~~ (covA_P45)*A_P5_T1 
    
    E_P1_T1 ~~ (covE_P12)*E_P2_T1 + (covE_P13)*E_P3_T1 + (covE_P14)*E_P4_T1 + (covE_P15)*E_P5_T1
    E_P2_T1 ~~ (covE_P23)*E_P3_T1 + (covE_P24)*E_P4_T1 + (covE_P25)*E_P5_T1
    E_P3_T1 ~~ (covE_P34)*E_P4_T1 + (covE_P35)*E_P5_T1
    E_P4_T1 ~~ (covE_P45)*E_P5_T1 
    
    A_P1_T2 ~~ (covA_P12)*A_P2_T2 + (covA_P13)*A_P3_T2 + (covA_P14)*A_P4_T2 + (covA_P15)*A_P5_T2
    A_P2_T2 ~~ (covA_P23)*A_P3_T2 + (covA_P24)*A_P4_T2 + (covA_P25)*A_P5_T2
    A_P3_T2 ~~ (covA_P34)*A_P4_T2 + (covA_P35)*A_P5_T2
    A_P4_T2 ~~ (covA_P45)*A_P5_T2 
    
    E_P1_T2 ~~ (covE_P12)*E_P2_T2 + (covE_P13)*E_P3_T2 + (covE_P14)*E_P4_T2 + (covE_P15)*E_P5_T2
    E_P2_T2 ~~ (covE_P23)*E_P3_T2 + (covE_P24)*E_P4_T2 + (covE_P25)*E_P5_T2
    E_P3_T2 ~~ (covE_P34)*E_P4_T2 + (covE_P35)*E_P5_T2
    E_P4_T2 ~~ (covE_P45)*E_P5_T2 
 
    #constrains
    cov_MZ_P12 == covA_P12
    cov_MZ_P13 == covA_P13
    cov_MZ_P14 == covA_P14
    cov_MZ_P15 == covA_P15
    cov_MZ_P23 == covA_P23
    cov_MZ_P24 == covA_P24
    cov_MZ_P25 == covA_P25
    cov_MZ_P34 == covA_P34
    cov_MZ_P35 == covA_P35
    cov_MZ_P45 == covA_P45
    
    cov_DZ_P12 == .5*covA_P12
    cov_DZ_P13 == .5*covA_P13
    cov_DZ_P14 == .5*covA_P14
    cov_DZ_P15 == .5*covA_P15
    cov_DZ_P23 == .5*covA_P23
    cov_DZ_P24 == .5*covA_P24
    cov_DZ_P25 == .5*covA_P25
    cov_DZ_P34 == .5*covA_P34
    cov_DZ_P35 == .5*covA_P35
    cov_DZ_P45 == .5*covA_P45
    
    #covariances cross-twin cross-trait
    A_P1_T1 ~~ c((cov_MZ_P12),(cov_MZ_P12),(cov_DZ_P12),(cov_DZ_P12),(cov_DZ_P12))*A_P2_T2 + c((cov_MZ_P13),(cov_MZ_P13),(cov_DZ_P13),(cov_DZ_P13),(cov_DZ_P13))*A_P3_T2 + c((cov_MZ_P14),(cov_MZ_P14),(cov_DZ_P14),(cov_DZ_P14),(cov_DZ_P14))*A_P4_T2 + c((cov_MZ_P15),(cov_MZ_P15),(cov_DZ_P15),(cov_DZ_P15),(cov_DZ_P15))*A_P5_T2
    A_P2_T1 ~~ c((cov_MZ_P12),(cov_MZ_P12),(cov_DZ_P12),(cov_DZ_P12),(cov_DZ_P12))*A_P1_T2 + c((cov_MZ_P23),(cov_MZ_P23),(cov_DZ_P23),(cov_DZ_P23),(cov_DZ_P23))*A_P3_T2 + c((cov_MZ_P24),(cov_MZ_P24),(cov_DZ_P24),(cov_DZ_P24),(cov_DZ_P24))*A_P4_T2 + c((cov_MZ_P25),(cov_MZ_P25),(cov_DZ_P25),(cov_DZ_P25),(cov_DZ_P25))*A_P5_T2
    A_P3_T1 ~~ c((cov_MZ_P13),(cov_MZ_P13),(cov_DZ_P13),(cov_DZ_P13),(cov_DZ_P13))*A_P1_T2 + c((cov_MZ_P23),(cov_MZ_P23),(cov_DZ_P23),(cov_DZ_P23),(cov_DZ_P23))*A_P2_T2 + c((cov_MZ_P34),(cov_MZ_P34),(cov_DZ_P34),(cov_DZ_P34),(cov_DZ_P34))*A_P4_T2 + c((cov_MZ_P35),(cov_MZ_P35),(cov_DZ_P35),(cov_DZ_P35),(cov_DZ_P35))*A_P5_T2
    A_P4_T1 ~~ c((cov_MZ_P14),(cov_MZ_P14),(cov_DZ_P14),(cov_DZ_P14),(cov_DZ_P14))*A_P1_T2 + c((cov_MZ_P24),(cov_MZ_P24),(cov_DZ_P24),(cov_DZ_P24),(cov_DZ_P24))*A_P2_T2 + c((cov_MZ_P34),(cov_MZ_P34),(cov_DZ_P34),(cov_DZ_P34),(cov_DZ_P34))*A_P3_T2 + c((cov_MZ_P45),(cov_MZ_P45),(cov_DZ_P45),(cov_DZ_P45),(cov_DZ_P45))*A_P5_T2
    A_P5_T1 ~~ c((cov_MZ_P15),(cov_MZ_P15),(cov_DZ_P15),(cov_DZ_P15),(cov_DZ_P15))*A_P1_T2 + c((cov_MZ_P25),(cov_MZ_P25),(cov_DZ_P25),(cov_DZ_P25),(cov_DZ_P25))*A_P2_T2 + c((cov_MZ_P35),(cov_MZ_P35),(cov_DZ_P35),(cov_DZ_P35),(cov_DZ_P35))*A_P3_T2 + c((cov_MZ_P45),(cov_MZ_P45),(cov_DZ_P45),(cov_DZ_P45),(cov_DZ_P45))*A_P4_T2
    
    
    E_P1_T1 ~~ 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P2_T1 ~~ 0*E_P1_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P3_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P4_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P5_T2
    E_P5_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2
    
    #set every other covariance to zero
    #within-trait A*G  
    A_P1_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2
    A_P1_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2
    
    A_P2_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2
    A_P2_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2
        
    A_P3_T1 ~~ 0*E_P3_T1 + 0*E_P3_T2
    A_P3_T2 ~~ 0*E_P3_T1 + 0*E_P3_T2

    A_P4_T1 ~~ 0*E_P4_T1 + 0*E_P4_T2
    A_P4_T2 ~~ 0*E_P4_T1 + 0*E_P4_T2
    
    A_P5_T1 ~~ 0*E_P5_T1 + 0*E_P5_T2
    A_P5_T2 ~~ 0*E_P5_T1 + 0*E_P5_T2
    

    #cross-trait A*G
    A_P1_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P1_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 

    A_P2_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P2_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    A_P3_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P3_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    A_P4_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P4_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    A_P5_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 
    A_P5_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 
    
    #estimates
    #phenotypic variance
    varP1 := (varA_P1+varE_P1)
    varP2 := (varA_P2+varE_P2)
    varP3 := (varA_P3+varE_P3)
    varP4 := (varA_P4+varE_P4)
    varP5 := (varA_P5+varE_P5)
    
    #compute proportion of covariance explained by genetics (bivariate heritabilities)
    bh2_P12 := covA_P12/ (covA_P12 + covE_P12)
    bh2_P13 := covA_P13/ (covA_P13 + covE_P13)
    bh2_P14 := covA_P14/ (covA_P14 + covE_P14)
    bh2_P15 := covA_P15/ (covA_P15 + covE_P15)
    bh2_P23 := covA_P23/ (covA_P23 + covE_P23)
    bh2_P24 := covA_P24/ (covA_P24 + covE_P24)
    bh2_P25 := covA_P25/ (covA_P25 + covE_P25)
    bh2_P34 := covA_P34/ (covA_P34 + covE_P34)
    bh2_P35 := covA_P35/ (covA_P35 + covE_P35)
    bh2_P45 := covA_P45/ (covA_P45 + covE_P45)
    
    #compute proportion of covariance explained by environement
    be2_P12 := covE_P12/ (covA_P12 + covE_P12)
    be2_P13 := covE_P13/ (covA_P13 + covE_P13)
    be2_P14 := covE_P14/ (covA_P14 + covE_P14)
    be2_P15 := covE_P15/ (covA_P15 + covE_P15)
    be2_P23 := covE_P23/ (covA_P23 + covE_P23)
    be2_P24 := covE_P24/ (covA_P24 + covE_P24)
    be2_P25 := covE_P25/ (covA_P25 + covE_P25)
    be2_P34 := covE_P34/ (covA_P34 + covE_P34)
    be2_P35 := covE_P35/ (covA_P35 + covE_P35)
    be2_P45 := covE_P45/ (covA_P45 + covE_P45)
    
    #compute genetic correlations
    sdA_P1 := sqrt(varA_P1)
    sdA_P2 := sqrt(varA_P2)
    sdA_P3 := sqrt(varA_P3)
    sdA_P4 := sqrt(varA_P4)
    sdA_P5 := sqrt(varA_P5)
    
    rA_P12 := covA_P12/(sdA_P1*sdA_P2)
    rA_P13 := covA_P13/(sdA_P1*sdA_P3)
    rA_P14 := covA_P14/(sdA_P1*sdA_P4)
    rA_P15 := covA_P15/(sdA_P1*sdA_P5)
    rA_P23 := covA_P23/(sdA_P2*sdA_P3)
    rA_P24 := covA_P24/(sdA_P2*sdA_P4)
    rA_P25 := covA_P25/(sdA_P2*sdA_P5)
    rA_P34 := covA_P34/(sdA_P3*sdA_P4)
    rA_P35 := covA_P35/(sdA_P3*sdA_P5)
    rA_P45 := covA_P45/(sdA_P4*sdA_P4)
    
    #compute environmental correlations
    sdE_P1 := sqrt(varE_P1)
    sdE_P2 := sqrt(varE_P2)
    sdE_P3 := sqrt(varE_P3)
    sdE_P4 := sqrt(varE_P4)
    sdE_P5 := sqrt(varE_P5)
    
    rE_P12 := covE_P12/(sdE_P1*sdE_P2)
    rE_P13 := covE_P13/(sdE_P1*sdE_P3)
    rE_P14 := covE_P14/(sdE_P1*sdE_P4)
    rE_P15 := covE_P15/(sdE_P1*sdE_P5)
    rE_P23 := covE_P23/(sdE_P2*sdE_P3)
    rE_P24 := covE_P24/(sdE_P2*sdE_P4)
    rE_P25 := covE_P25/(sdE_P2*sdE_P5)
    rE_P34 := covE_P34/(sdE_P3*sdE_P4)
    rE_P35 := covE_P35/(sdE_P3*sdE_P5)
    rE_P45 := covE_P45/(sdE_P4*sdE_P5)
    
    #twin-h2
    h2_P1 := varA_P1/ varP1
    h2_P2 := varA_P2/ varP2
    h2_P3 := varA_P3/ varP3
    h2_P4 := varA_P4/ varP4
    h2_P5 := varA_P5/ varP5
    
    #rA difference
    d_P123 := covA_P12 - covA_P13
    d_P124 := covA_P12 - covA_P14
    d_P125 := covA_P12 - covA_P15
    d_P134 := covA_P13 - covA_P14
    d_P135 := covA_P13 - covA_P15
    d_P145 := covA_P14 - covA_P15
    
  "


hiapm_AE_mod <-"
    #means
    P1_T1~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_f_P1)*1 + b_P1*age
    P1_T2~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_m_P1)*1 + b_P1*age
    P2_T1~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_f_P2)*1 + b_P2*age
    P2_T2~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_m_P2)*1 + b_P2*age
    P3_T1~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_f_P3)*1 + b_P3*age
    P3_T2~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_m_P3)*1 + b_P3*age
    P4_T1~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4, mean_f_P4)*1 + b_P4*age
    P4_T2~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4, mean_m_P4)*1 + b_P4*age
    P5_T1~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5, mean_f_P5)*1 + b_P5*age
    P5_T2~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5, mean_m_P5)*1 + b_P5*age
    
    #constrain covariate variances to be equal across groups
    age ~~ var_age*age
    
    #shared A components
    A_1=~ NA*P1_T1 + A_P1*P1_T1 + A_P2*P2_T1 + A_P3*P3_T1 + A_P4*P4_T1 + A_P5*P5_T1
    A_2=~ NA*P1_T2 + A_P1*P1_T2 + A_P2*P2_T2 + A_P3*P3_T2 + A_P4*P4_T2 + A_P5*P5_T2
    
    #variances
    A_1 ~~ 1*A_1
    A_2 ~~ 1*A_2
    
    #AE components (unique)
    A_P1_T1=~ P1_T1
    A_P1_T2=~ P1_T2
    A_P2_T1=~ P2_T1
    A_P2_T2=~ P2_T2
    A_P3_T1=~ P3_T1
    A_P3_T2=~ P3_T2   
    A_P4_T1=~ P4_T1
    A_P4_T2=~ P4_T2
    A_P5_T1=~ P5_T1
    A_P5_T2=~ P5_T2
    
    E_P1_T1=~ P1_T1
    E_P1_T2=~ P1_T2
    E_P2_T1=~ P2_T1
    E_P2_T2=~ P2_T2
    E_P3_T1=~ P3_T1
    E_P3_T2=~ P3_T2   
    E_P4_T1=~ P4_T1
    E_P4_T2=~ P4_T2
    E_P5_T1=~ P5_T1
    E_P5_T2=~ P5_T2

    #variances
    A_P1_T1 ~~ varA_P1*A_P1_T1
    A_P1_T2 ~~ varA_P1*A_P1_T2
    A_P2_T1 ~~ varA_P2*A_P2_T1
    A_P2_T2 ~~ varA_P2*A_P2_T2
    A_P3_T1 ~~ varA_P3*A_P3_T1
    A_P3_T2 ~~ varA_P3*A_P3_T2
    A_P4_T1 ~~ varA_P4*A_P4_T1
    A_P4_T2 ~~ varA_P4*A_P4_T2
    A_P5_T1 ~~ varA_P5*A_P5_T1
    A_P5_T2 ~~ varA_P5*A_P5_T2

    E_P1_T1 ~~ varE_P1*E_P1_T1
    E_P1_T2 ~~ varE_P1*E_P1_T2
    E_P2_T1 ~~ varE_P2*E_P2_T1
    E_P2_T2 ~~ varE_P2*E_P2_T2
    E_P3_T1 ~~ varE_P3*E_P3_T1
    E_P3_T2 ~~ varE_P3*E_P3_T2
    E_P4_T1 ~~ varE_P4*E_P4_T1
    E_P4_T2 ~~ varE_P4*E_P4_T2
    E_P5_T1 ~~ varE_P5*E_P5_T1
    E_P5_T2 ~~ varE_P5*E_P5_T2
    
    #fix remaining variance to 0
    P1_T1~~0*P1_T1
    P1_T2~~0*P1_T2
    P2_T1~~0*P2_T1
    P2_T2~~0*P2_T2
    P3_T1~~0*P3_T1
    P3_T2~~0*P3_T2
    P4_T1~~0*P4_T1
    P4_T2~~0*P4_T2
    P5_T1~~0*P5_T1
    P5_T2~~0*P5_T2

    #constraints (twin pair covariances, p = phenotype)
    cov_MZ_A == 1
    cov_DZ_A == .5
    
    cov_MZ_P1 == varA_P1
    cov_MZ_P2 == varA_P2
    cov_MZ_P3 == varA_P3
    cov_MZ_P4 == varA_P4
    cov_MZ_P5 == varA_P5

    cov_DZ_P1 == .5*varA_P1
    cov_DZ_P2 == .5*varA_P2
    cov_DZ_P3 == .5*varA_P3
    cov_DZ_P4 == .5*varA_P4
    cov_DZ_P5 == .5*varA_P5
    
    #covariances (shared)
    A_1 ~~ c((cov_MZ_A),(cov_MZ_A),(cov_DZ_A),(cov_DZ_A),(cov_DZ_A))*A_2
    
    #covariances cross-trait within twin (unique)
    A_P1_T1 ~~ c((cov_MZ_P1),(cov_MZ_P1),(cov_DZ_P1),(cov_DZ_P1),(cov_DZ_P1))*A_P1_T2
    A_P2_T1 ~~ c((cov_MZ_P2),(cov_MZ_P2),(cov_DZ_P2),(cov_DZ_P2),(cov_DZ_P2))*A_P2_T2
    A_P3_T1 ~~ c((cov_MZ_P3),(cov_MZ_P3),(cov_DZ_P3),(cov_DZ_P3),(cov_DZ_P3))*A_P3_T2 
    A_P4_T1 ~~ c((cov_MZ_P4),(cov_MZ_P4),(cov_DZ_P4),(cov_DZ_P4),(cov_DZ_P4))*A_P4_T2
    A_P5_T1 ~~ c((cov_MZ_P5),(cov_MZ_P5),(cov_DZ_P5),(cov_DZ_P5),(cov_DZ_P5))*A_P5_T2
    

    #covariances within-twin cross-trait
    E_P1_T1 ~~ (covE_P12)*E_P2_T1 + (covE_P13)*E_P3_T1 + (covE_P14)*E_P4_T1 + (covE_P15)*E_P5_T1
    E_P2_T1 ~~ (covE_P23)*E_P3_T1 + (covE_P24)*E_P4_T1 + (covE_P25)*E_P5_T1
    E_P3_T1 ~~ (covE_P34)*E_P4_T1 + (covE_P35)*E_P5_T1
    E_P4_T1 ~~ (covE_P45)*E_P5_T1 
    
    E_P1_T2 ~~ (covE_P12)*E_P2_T2 + (covE_P13)*E_P3_T2 + (covE_P14)*E_P4_T2 + (covE_P15)*E_P5_T2
    E_P2_T2 ~~ (covE_P23)*E_P3_T2 + (covE_P24)*E_P4_T2 + (covE_P25)*E_P5_T2
    E_P3_T2 ~~ (covE_P34)*E_P4_T2 + (covE_P35)*E_P5_T2
    E_P4_T2 ~~ (covE_P45)*E_P5_T2 
 
    #cross-trait A*G (unique)
    A_P1_T1 ~~ 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P1_T2 ~~ 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P3_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P3_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P2_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P2_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P4_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P4_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P5_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2
    A_P5_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2
   
    #cross-trait A (unique)
    A_P1_T1 ~~ 0*A_P2_T2 + 0*A_P3_T2 + 0*A_P4_T2 + 0*A_P5_T2 + 0*A_P2_T1 + 0*A_P3_T1 + 0*A_P4_T1 + 0*A_P5_T1
    A_P2_T1 ~~ 0*A_P1_T2 + 0*A_P3_T2 + 0*A_P4_T2 + 0*A_P5_T2 + 0*A_P3_T1 + 0*A_P4_T1 + 0*A_P5_T1
    A_P3_T1 ~~ 0*A_P1_T2 + 0*A_P2_T2 + 0*A_P4_T2 + 0*A_P5_T2 + 0*A_P4_T1 + 0*A_P5_T1
    A_P4_T1 ~~ 0*A_P1_T2 + 0*A_P2_T2 + 0*A_P3_T2 + 0*A_P5_T2 + 0*A_P5_T1
    A_P5_T1 ~~ 0*A_P1_T2 + 0*A_P2_T2 + 0*A_P3_T2 + 0*A_P4_T2 
    
    A_P1_T2 ~~ 0*A_P2_T2 + 0*A_P3_T2 + 0*A_P4_T2 + 0*A_P5_T2 
    A_P2_T2 ~~ 0*A_P3_T2 + 0*A_P4_T2 + 0*A_P5_T2
    A_P3_T2 ~~ 0*A_P4_T2 + 0*A_P5_T2
    A_P4_T2 ~~ 0*A_P5_T2
  
    
    #cross trait E (unique)
    E_P1_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2 
    E_P2_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P3_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2 
    E_P4_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2 
    E_P5_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2 
    
    #unique*shared
    A_1 ~~ 0*A_P1_T1 + 0*A_P1_T2 + 0*A_P2_T1 + 0*A_P2_T2 + 0*A_P3_T1 + 0*A_P3_T2 + 0*A_P4_T1 + 0*A_P4_T2 + 0*A_P5_T1 + 0*A_P5_T2  + 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_2 ~~ 0*A_P1_T1 + 0*A_P1_T2 + 0*A_P2_T1 + 0*A_P2_T2 + 0*A_P3_T1 + 0*A_P3_T2 + 0*A_P4_T1 + 0*A_P4_T2 + 0*A_P5_T1 + 0*A_P5_T2  + 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
     #set every other covariance to zero
    #within-trait A*G  
    A_P1_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2
    A_P1_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2
    A_P2_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2
    A_P2_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2
    A_P3_T1 ~~ 0*E_P3_T1 + 0*E_P3_T2
    A_P3_T2 ~~ 0*E_P3_T1 + 0*E_P3_T2
    A_P4_T1 ~~ 0*E_P4_T1 + 0*E_P4_T2
    A_P4_T2 ~~ 0*E_P4_T1 + 0*E_P4_T2
    A_P5_T1 ~~ 0*E_P5_T1 + 0*E_P5_T2
    A_P5_T2 ~~ 0*E_P5_T1 + 0*E_P5_T2
    
    
    #phenotypic variance
    #variances
    varP1 := (A_P1^2) + varA_P1 + varE_P1
    varP2 := (A_P2^2) + varA_P2 + varE_P2
    varP3 := (A_P3^2) + varA_P3 + varE_P3
    varP4 := (A_P4^2) + varA_P4 + varE_P4
    varP5 := (A_P5^2) + varA_P5 + varE_P5

    #compute environmental correlations
    sdE_P1 := sqrt(varE_P1)
    sdE_P2 := sqrt(varE_P2)
    sdE_P3 := sqrt(varE_P3)
    sdE_P4 := sqrt(varE_P4)
    sdE_P5 := sqrt(varE_P5)

    rE_P12 := covE_P12/(sdE_P1*sdE_P2)
    rE_P13 := covE_P13/(sdE_P1*sdE_P3)
    rE_P14 := covE_P14/(sdE_P1*sdE_P4)
    rE_P15 := covE_P15/(sdE_P1*sdE_P5)
    rE_P23 := covE_P23/(sdE_P2*sdE_P3)
    rE_P24 := covE_P24/(sdE_P2*sdE_P4)
    rE_P25 := covE_P25/(sdE_P2*sdE_P5)
    rE_P34 := covE_P34/(sdE_P3*sdE_P4)
    rE_P35 := covE_P35/(sdE_P3*sdE_P5)
    rE_P45 := covE_P45/(sdE_P4*sdE_P5)

    #shared
    A_P1_shared := (A_P1^2)/varP1
    A_P2_shared := (A_P2^2)/varP2
    A_P3_shared := (A_P3^2)/varP3
    A_P4_shared := (A_P4^2)/varP4
    A_P5_shared := (A_P5^2)/varP5
    
    #unique
    A_P1_unique := varA_P1/varP1
    A_P2_unique := varA_P2/varP2
    A_P3_unique := varA_P3/varP3
    A_P4_unique := varA_P4/varP4
    A_P5_unique := varA_P5/varP5
    
    #twin-h2
    h2p1 := ((A_P1^2) + varA_P1)/ varP1
    h2p2 := ((A_P2^2) + varA_P2)/ varP2
    h2p3 := ((A_P3^2) + varA_P3)/ varP3
    h2p4 := ((A_P4^2) + varA_P4)/ varP4
    h2p5 := ((A_P5^2) + varA_P5)/ varP5
  "

hiepm_AE_mod <-"
    #means
    P1_T1~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_f_P1)*1 + b_P1*age
    P1_T2~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_m_P1)*1 + b_P1*age
    P2_T1~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_f_P2)*1 + b_P2*age
    P2_T2~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_m_P2)*1 + b_P2*age
    P3_T1~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_f_P3)*1 + b_P3*age
    P3_T2~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_m_P3)*1 + b_P3*age
    P4_T1~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4, mean_f_P4)*1 + b_P4*age
    P4_T2~c(mean_f_P4, mean_m_P4, mean_f_P4, mean_m_P4, mean_m_P4)*1 + b_P4*age
    P5_T1~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5, mean_f_P5)*1 + b_P5*age
    P5_T2~c(mean_f_P5, mean_m_P5, mean_f_P5, mean_m_P5, mean_m_P5)*1 + b_P5*age
    
    #constrain covariate variances to be equal across groups
    age ~~ var_age*age
    
    #shared E components
    E_1=~ NA*P1_T1 + E_P1*P1_T1 + E_P2*P2_T1 + E_P3*P3_T1 + E_P4*P4_T1 + E_P5*P5_T1
    E_2=~ NA*P1_T2 + E_P1*P1_T2 + E_P2*P2_T2 + E_P3*P3_T2 + E_P4*P4_T2 + E_P5*P5_T2
    
    #variances
    E_1 ~~ 1*E_1
    E_2 ~~ 1*E_2
    
    #AE components
    A_P1_T1=~ P1_T1
    A_P1_T2=~ P1_T2
    A_P2_T1=~ P2_T1
    A_P2_T2=~ P2_T2
    A_P3_T1=~ P3_T1
    A_P3_T2=~ P3_T2   
    A_P4_T1=~ P4_T1
    A_P4_T2=~ P4_T2
    A_P5_T1=~ P5_T1
    A_P5_T2=~ P5_T2

    E_P1_T1=~ P1_T1
    E_P1_T2=~ P1_T2
    E_P2_T1=~ P2_T1
    E_P2_T2=~ P2_T2
    E_P3_T1=~ P3_T1
    E_P3_T2=~ P3_T2   
    E_P4_T1=~ P4_T1
    E_P4_T2=~ P4_T2
    E_P5_T1=~ P5_T1
    E_P5_T2=~ P5_T2

    #variances
    A_P1_T1 ~~ varA_P1*A_P1_T1
    A_P1_T2 ~~ varA_P1*A_P1_T2
    A_P2_T1 ~~ varA_P2*A_P2_T1
    A_P2_T2 ~~ varA_P2*A_P2_T2
    A_P3_T1 ~~ varA_P3*A_P3_T1
    A_P3_T2 ~~ varA_P3*A_P3_T2
    A_P4_T1 ~~ varA_P4*A_P4_T1
    A_P4_T2 ~~ varA_P4*A_P4_T2
    A_P5_T1 ~~ varA_P5*A_P5_T1
    A_P5_T2 ~~ varA_P5*A_P5_T2

    E_P1_T1 ~~ varE_P1*E_P1_T1
    E_P1_T2 ~~ varE_P1*E_P1_T2
    E_P2_T1 ~~ varE_P2*E_P2_T1
    E_P2_T2 ~~ varE_P2*E_P2_T2
    E_P3_T1 ~~ varE_P3*E_P3_T1
    E_P3_T2 ~~ varE_P3*E_P3_T2
    E_P4_T1 ~~ varE_P4*E_P4_T1
    E_P4_T2 ~~ varE_P4*E_P4_T2
    E_P5_T1 ~~ varE_P5*E_P5_T1
    E_P5_T2 ~~ varE_P5*E_P5_T2

    #fix remaining variance to 0
    P1_T1~~0*P1_T1
    P1_T2~~0*P1_T2
    P2_T1~~0*P2_T1
    P2_T2~~0*P2_T2
    P3_T1~~0*P3_T1
    P3_T2~~0*P3_T2
    P4_T1~~0*P4_T1
    P4_T2~~0*P4_T2
    P5_T1~~0*P5_T1
    P5_T2~~0*P5_T2
    
    
    #constraints (twin pair covariances, p = phenotype)
    cov_MZ_P1 == varA_P1
    cov_MZ_P2 == varA_P2
    cov_MZ_P3 == varA_P3
    cov_MZ_P4 == varA_P4
    cov_MZ_P5 == varA_P5

    cov_DZ_P1 == .5*varA_P1
    cov_DZ_P2 == .5*varA_P2
    cov_DZ_P3 == .5*varA_P3
    cov_DZ_P4 == .5*varA_P4
    cov_DZ_P5 == .5*varA_P5
    
    #covariances cross-twin wihtin-trait
    A_P1_T1 ~~ c((cov_MZ_P1),(cov_MZ_P1),(cov_DZ_P1),(cov_DZ_P1),(cov_DZ_P1))*A_P1_T2
    A_P2_T1 ~~ c((cov_MZ_P2),(cov_MZ_P2),(cov_DZ_P2),(cov_DZ_P2),(cov_DZ_P2))*A_P2_T2
    A_P3_T1 ~~ c((cov_MZ_P3),(cov_MZ_P3),(cov_DZ_P3),(cov_DZ_P3),(cov_DZ_P3))*A_P3_T2 
    A_P4_T1 ~~ c((cov_MZ_P4),(cov_MZ_P4),(cov_DZ_P4),(cov_DZ_P4),(cov_DZ_P4))*A_P4_T2
    A_P5_T1 ~~ c((cov_MZ_P5),(cov_MZ_P5),(cov_DZ_P5),(cov_DZ_P5),(cov_DZ_P5))*A_P5_T2
    
    E_P1_T1 ~~ 0*E_P1_T2
    E_P2_T1 ~~ 0*E_P2_T2
    E_P3_T1 ~~ 0*E_P3_T2
    E_P4_T1 ~~ 0*E_P4_T2
    E_P5_T1 ~~ 0*E_P5_T2

    #covariances within-twin cross-trait
    A_P1_T1 ~~ (covA_P12)*A_P2_T1 + (covA_P13)*A_P3_T1 + (covA_P14)*A_P4_T1 + (covA_P15)*A_P5_T1
    A_P2_T1 ~~ (covA_P23)*A_P3_T1 + (covA_P24)*A_P4_T1 + (covA_P25)*A_P5_T1
    A_P3_T1 ~~ (covA_P34)*A_P4_T1 + (covA_P35)*A_P5_T1
    A_P4_T1 ~~ (covA_P45)*A_P5_T1 
    
    E_P1_T1 ~~ 0*E_P2_T1 + 0*E_P3_T1 + 0*E_P4_T1 + 0*E_P5_T1
    E_P2_T1 ~~ 0*E_P3_T1 + 0*E_P4_T1 + 0*E_P5_T1
    E_P3_T1 ~~ 0*E_P4_T1 + 0*E_P5_T1
    E_P4_T1 ~~ 0*E_P5_T1 
    
    A_P1_T2 ~~ (covA_P12)*A_P2_T2 + (covA_P13)*A_P3_T2 + (covA_P14)*A_P4_T2 + (covA_P15)*A_P5_T2
    A_P2_T2 ~~ (covA_P23)*A_P3_T2 + (covA_P24)*A_P4_T2 + (covA_P25)*A_P5_T2
    A_P3_T2 ~~ (covA_P34)*A_P4_T2 + (covA_P35)*A_P5_T2
    A_P4_T2 ~~ (covA_P45)*A_P5_T2 
    
    E_P1_T2 ~~ 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P2_T2 ~~ 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P3_T2 ~~ 0*E_P4_T2 + 0*E_P5_T2
    E_P4_T2 ~~ 0*E_P5_T2 
 
    #constrains
    cov_MZ_P12 == covA_P12
    cov_MZ_P13 == covA_P13
    cov_MZ_P14 == covA_P14
    cov_MZ_P15 == covA_P15
    cov_MZ_P23 == covA_P23
    cov_MZ_P24 == covA_P24
    cov_MZ_P25 == covA_P25
    cov_MZ_P34 == covA_P34
    cov_MZ_P35 == covA_P35
    cov_MZ_P45 == covA_P45
    
    cov_DZ_P12 == .5*covA_P12
    cov_DZ_P13 == .5*covA_P13
    cov_DZ_P14 == .5*covA_P14
    cov_DZ_P15 == .5*covA_P15
    cov_DZ_P23 == .5*covA_P23
    cov_DZ_P24 == .5*covA_P24
    cov_DZ_P25 == .5*covA_P25
    cov_DZ_P34 == .5*covA_P34
    cov_DZ_P35 == .5*covA_P35
    cov_DZ_P45 == .5*covA_P45
    
    #covariances cross-twin cross-trait
    A_P1_T1 ~~ c((cov_MZ_P12),(cov_MZ_P12),(cov_DZ_P12),(cov_DZ_P12),(cov_DZ_P12))*A_P2_T2 + c((cov_MZ_P13),(cov_MZ_P13),(cov_DZ_P13),(cov_DZ_P13),(cov_DZ_P13))*A_P3_T2 + c((cov_MZ_P14),(cov_MZ_P14),(cov_DZ_P14),(cov_DZ_P14),(cov_DZ_P14))*A_P4_T2 + c((cov_MZ_P15),(cov_MZ_P15),(cov_DZ_P15),(cov_DZ_P15),(cov_DZ_P15))*A_P5_T2
    A_P2_T1 ~~ c((cov_MZ_P12),(cov_MZ_P12),(cov_DZ_P12),(cov_DZ_P12),(cov_DZ_P12))*A_P1_T2 + c((cov_MZ_P23),(cov_MZ_P23),(cov_DZ_P23),(cov_DZ_P23),(cov_DZ_P23))*A_P3_T2 + c((cov_MZ_P24),(cov_MZ_P24),(cov_DZ_P24),(cov_DZ_P24),(cov_DZ_P24))*A_P4_T2 + c((cov_MZ_P25),(cov_MZ_P25),(cov_DZ_P25),(cov_DZ_P25),(cov_DZ_P25))*A_P5_T2
    A_P3_T1 ~~ c((cov_MZ_P13),(cov_MZ_P13),(cov_DZ_P13),(cov_DZ_P13),(cov_DZ_P13))*A_P1_T2 + c((cov_MZ_P23),(cov_MZ_P23),(cov_DZ_P23),(cov_DZ_P23),(cov_DZ_P23))*A_P2_T2 + c((cov_MZ_P34),(cov_MZ_P34),(cov_DZ_P34),(cov_DZ_P34),(cov_DZ_P34))*A_P4_T2 + c((cov_MZ_P35),(cov_MZ_P35),(cov_DZ_P35),(cov_DZ_P35),(cov_DZ_P35))*A_P5_T2
    A_P4_T1 ~~ c((cov_MZ_P14),(cov_MZ_P14),(cov_DZ_P14),(cov_DZ_P14),(cov_DZ_P14))*A_P1_T2 + c((cov_MZ_P24),(cov_MZ_P24),(cov_DZ_P24),(cov_DZ_P24),(cov_DZ_P24))*A_P2_T2 + c((cov_MZ_P34),(cov_MZ_P34),(cov_DZ_P34),(cov_DZ_P34),(cov_DZ_P34))*A_P3_T2 + c((cov_MZ_P45),(cov_MZ_P45),(cov_DZ_P45),(cov_DZ_P45),(cov_DZ_P45))*A_P5_T2
    A_P5_T1 ~~ c((cov_MZ_P15),(cov_MZ_P15),(cov_DZ_P15),(cov_DZ_P15),(cov_DZ_P15))*A_P1_T2 + c((cov_MZ_P25),(cov_MZ_P25),(cov_DZ_P25),(cov_DZ_P25),(cov_DZ_P25))*A_P2_T2 + c((cov_MZ_P35),(cov_MZ_P35),(cov_DZ_P35),(cov_DZ_P35),(cov_DZ_P35))*A_P3_T2 + c((cov_MZ_P45),(cov_MZ_P45),(cov_DZ_P45),(cov_DZ_P45),(cov_DZ_P45))*A_P4_T2
    
    E_P1_T1 ~~ 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P2_T1 ~~ 0*E_P1_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P3_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P4_T2 + 0*E_P5_T2
    E_P4_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P5_T2
    E_P5_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2
    
    #set every other covariance to zero
    #within-trait A*G  
    A_P1_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2
    A_P1_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2
    
    A_P2_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2
    A_P2_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2
        
    A_P3_T1 ~~ 0*E_P3_T1 + 0*E_P3_T2
    A_P3_T2 ~~ 0*E_P3_T1 + 0*E_P3_T2

    A_P4_T1 ~~ 0*E_P4_T1 + 0*E_P4_T2
    A_P4_T2 ~~ 0*E_P4_T1 + 0*E_P4_T2
    
    A_P5_T1 ~~ 0*E_P5_T1 + 0*E_P5_T2
    A_P5_T2 ~~ 0*E_P5_T1 + 0*E_P5_T2
    

    #cross-trait A*G
    A_P1_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P1_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 

    A_P2_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2
    A_P2_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    A_P3_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P3_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    A_P4_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    A_P4_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    A_P5_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 
    A_P5_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 
    
    #unique*shared
    E_1 ~~ 0*A_P1_T1 + 0*A_P1_T2 + 0*A_P2_T1 + 0*A_P2_T2 + 0*A_P3_T1 + 0*A_P3_T2 + 0*A_P4_T1 + 0*A_P4_T2 + 0*A_P5_T1 + 0*A_P5_T2  + 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    E_2 ~~ 0*A_P1_T1 + 0*A_P1_T2 + 0*A_P2_T1 + 0*A_P2_T2 + 0*A_P3_T1 + 0*A_P3_T2 + 0*A_P4_T1 + 0*A_P4_T2 + 0*A_P5_T1 + 0*A_P5_T2  + 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 

    #shared
    E_1 ~~ 0*E_2

    #phenotypic variance
    #variances
    varP1 := (E_P1^2) + varA_P1 + varE_P1
    varP2 := (E_P2^2) + varA_P2 + varE_P2
    varP3 := (E_P3^2) + varA_P3 + varE_P3
    varP4 := (E_P4^2) + varA_P4 + varE_P4
    varP5 := (E_P5^2) + varA_P5 + varE_P5

    #compute genetic correlations
    sdA_P1 := sqrt(varA_P1)
    sdA_P2 := sqrt(varA_P2)
    sdA_P3 := sqrt(varA_P3)
    sdA_P4 := sqrt(varA_P4)
    sdA_P5 := sqrt(varA_P5)

    rA_P12 := covA_P12/(sdA_P1*sdA_P2)
    rA_P13 := covA_P13/(sdA_P1*sdA_P3)
    rA_P14 := covA_P14/(sdA_P1*sdA_P4)
    rA_P15 := covA_P15/(sdA_P1*sdA_P5)
    rA_P23 := covA_P23/(sdA_P2*sdA_P3)
    rA_P24 := covA_P24/(sdA_P2*sdA_P4)
    rA_P25 := covA_P25/(sdA_P2*sdA_P5)
    rA_P34 := covA_P34/(sdA_P3*sdA_P4)
    rA_P35 := covA_P35/(sdA_P3*sdA_P5)
    rA_P45 := covA_P45/(sdA_P4*sdA_P5)

    #shared
    E_P1_shared := (E_P1^2)/varP1
    E_P2_shared := (E_P2^2)/varP2
    E_P3_shared := (E_P3^2)/varP3
    E_P4_shared := (E_P4^2)/varP4
    E_P5_shared := (E_P5^2)/varP5
    
    #unique
    E_P1_unique := varE_P1/varP1
    E_P2_unique := varE_P2/varP2
    E_P3_unique := varE_P3/varP3
    E_P4_unique := varE_P4/varP4
    E_P5_unique := varE_P5/varP5
    
    #twin-h2
    h2p1 := (varA_P1)/ varP1
    h2p2 := (varA_P2)/ varP2
    h2p3 := (varA_P3)/ varP3
    h2p4 := (varA_P4)/ varP4
    h2p5 := (varA_P5)/ varP5
    
  "