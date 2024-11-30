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
    P6_T1~c(mean_f_P6, mean_m_P6, mean_f_P6, mean_m_P6, mean_f_P6)*1 + b_P5*age
    P6_T2~c(mean_f_P6, mean_m_P6, mean_f_P6, mean_m_P6, mean_m_P6)*1 + b_P5*age
    
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
    A_P6_T1=~ P6_T1
    A_P6_T2=~ P6_T2

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
    E_P6_T1=~ P6_T1
    E_P6_T2=~ P6_T2

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
    A_P6_T1 ~~ varA_P6*A_P6_T1
    A_P6_T2 ~~ varA_P6*A_P6_T2

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
    E_P6_T1 ~~ varE_P6*E_P6_T1
    E_P6_T2 ~~ varE_P6*E_P6_T2

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
    P6_T1~~0*P6_T1
    P6_T2~~0*P6_T2
    
    
    #constraints (twin pair covariances, p = phenotype)
    cov_MZ_P1 == varA_P1
    cov_MZ_P2 == varA_P2
    cov_MZ_P3 == varA_P3
    cov_MZ_P4 == varA_P4
    cov_MZ_P5 == varA_P5
    cov_MZ_P6 == varA_P6

    cov_DZ_P1 == .5*varA_P1
    cov_DZ_P2 == .5*varA_P2
    cov_DZ_P3 == .5*varA_P3
    cov_DZ_P4 == .5*varA_P4
    cov_DZ_P5 == .5*varA_P5
    cov_DZ_P6 == .5*varA_P6
    
    #covariances cross-twin wihtin-trait
    A_P1_T1 ~~ c((cov_MZ_P1),(cov_MZ_P1),(cov_DZ_P1),(cov_DZ_P1),(cov_DZ_P1))*A_P1_T2
    A_P2_T1 ~~ c((cov_MZ_P2),(cov_MZ_P2),(cov_DZ_P2),(cov_DZ_P2),(cov_DZ_P2))*A_P2_T2
    A_P3_T1 ~~ c((cov_MZ_P3),(cov_MZ_P3),(cov_DZ_P3),(cov_DZ_P3),(cov_DZ_P3))*A_P3_T2 
    A_P4_T1 ~~ c((cov_MZ_P4),(cov_MZ_P4),(cov_DZ_P4),(cov_DZ_P4),(cov_DZ_P4))*A_P4_T2
    A_P5_T1 ~~ c((cov_MZ_P5),(cov_MZ_P5),(cov_DZ_P5),(cov_DZ_P5),(cov_DZ_P5))*A_P5_T2
    A_P6_T1 ~~ c((cov_MZ_P6),(cov_MZ_P6),(cov_DZ_P6),(cov_DZ_P6),(cov_DZ_P6))*A_P6_T2
    
    E_P1_T1 ~~ 0*E_P1_T2
    E_P2_T1 ~~ 0*E_P2_T2
    E_P3_T1 ~~ 0*E_P3_T2
    E_P4_T1 ~~ 0*E_P4_T2
    E_P5_T1 ~~ 0*E_P5_T2
    E_P6_T1 ~~ 0*E_P6_T2

    #covariances within-twin cross-trait
    A_P1_T1 ~~ (covA_P12)*A_P2_T1 + (covA_P13)*A_P3_T1 + (covA_P14)*A_P4_T1 + (covA_P15)*A_P5_T1 + (covA_P16)*A_P6_T1
    A_P2_T1 ~~ (covA_P23)*A_P3_T1 + (covA_P24)*A_P4_T1 + (covA_P25)*A_P5_T1 + (covA_P26)*A_P6_T1
    A_P3_T1 ~~ (covA_P34)*A_P4_T1 + (covA_P35)*A_P5_T1 + (covA_P36)*A_P6_T1
    A_P4_T1 ~~ (covA_P45)*A_P5_T1 + (covA_P46)*A_P6_T1
    A_P5_T1 ~~ (covA_P56)*A_P6_T1
    
    E_P1_T1 ~~ (covE_P12)*E_P2_T1 + (covE_P13)*E_P3_T1 + (covE_P14)*E_P4_T1 + (covE_P15)*E_P5_T1 + (covE_P16)*E_P6_T1
    E_P2_T1 ~~ (covE_P23)*E_P3_T1 + (covE_P24)*E_P4_T1 + (covE_P25)*E_P5_T1 + (covE_P26)*E_P6_T1
    E_P3_T1 ~~ (covE_P34)*E_P4_T1 + (covE_P35)*E_P5_T1 + (covE_P36)*E_P6_T1
    E_P4_T1 ~~ (covE_P45)*E_P5_T1 + (covE_P46)*E_P6_T1
    E_P5_T1 ~~ (covE_P56)*E_P6_T1
    
    A_P1_T2 ~~ (covA_P12)*A_P2_T2 + (covA_P13)*A_P3_T2 + (covA_P14)*A_P4_T2 + (covA_P15)*A_P5_T2 + (covA_P16)*A_P6_T2
    A_P2_T2 ~~ (covA_P23)*A_P3_T2 + (covA_P24)*A_P4_T2 + (covA_P25)*A_P5_T2 + (covA_P26)*A_P6_T2
    A_P3_T2 ~~ (covA_P34)*A_P4_T2 + (covA_P35)*A_P5_T2 + (covA_P36)*A_P6_T2
    A_P4_T2 ~~ (covA_P45)*A_P5_T2 + (covA_P46)*A_P6_T2
    A_P5_T2 ~~ (covA_P56)*A_P6_T2 
    
    E_P1_T2 ~~ (covE_P12)*E_P2_T2 + (covE_P13)*E_P3_T2 + (covE_P14)*E_P4_T2 + (covE_P15)*E_P5_T2 + (covE_P16)*E_P6_T2
    E_P2_T2 ~~ (covE_P23)*E_P3_T2 + (covE_P24)*E_P4_T2 + (covE_P25)*E_P5_T2 + (covE_P26)*E_P6_T2
    E_P3_T2 ~~ (covE_P34)*E_P4_T2 + (covE_P35)*E_P5_T2 + (covE_P36)*E_P6_T2
    E_P4_T2 ~~ (covE_P45)*E_P5_T2 + (covE_P46)*E_P6_T2
    E_P5_T2 ~~ (covE_P56)*E_P6_T2
 
    #constrains
    cov_MZ_P12 == covA_P12
    cov_MZ_P13 == covA_P13
    cov_MZ_P14 == covA_P14
    cov_MZ_P15 == covA_P15
    cov_MZ_P16 == covA_P16
    cov_MZ_P23 == covA_P23
    cov_MZ_P24 == covA_P24
    cov_MZ_P25 == covA_P25
    cov_MZ_P26 == covA_P26
    cov_MZ_P34 == covA_P34
    cov_MZ_P35 == covA_P35
    cov_MZ_P36 == covA_P36
    cov_MZ_P45 == covA_P45
    cov_MZ_P46 == covA_P46
    cov_MZ_P56 == covA_P56
    
    cov_DZ_P12 == .5*covA_P12
    cov_DZ_P13 == .5*covA_P13
    cov_DZ_P14 == .5*covA_P14
    cov_DZ_P15 == .5*covA_P15
    cov_DZ_P16 == .5*covA_P16
    cov_DZ_P23 == .5*covA_P23
    cov_DZ_P24 == .5*covA_P24
    cov_DZ_P25 == .5*covA_P25
    cov_DZ_P26 == .5*covA_P26
    cov_DZ_P34 == .5*covA_P34
    cov_DZ_P35 == .5*covA_P35
    cov_DZ_P36 == .5*covA_P36
    cov_DZ_P45 == .5*covA_P45
    cov_DZ_P46 == .5*covA_P46
    cov_DZ_P56 == .5*covA_P56
    
    #covariances cross-twin cross-trait
    A_P1_T1 ~~ c((cov_MZ_P12),(cov_MZ_P12),(cov_DZ_P12),(cov_DZ_P12),(cov_DZ_P12))*A_P2_T2 + c((cov_MZ_P13),(cov_MZ_P13),(cov_DZ_P13),(cov_DZ_P13),(cov_DZ_P13))*A_P3_T2 + c((cov_MZ_P14),(cov_MZ_P14),(cov_DZ_P14),(cov_DZ_P14),(cov_DZ_P14))*A_P4_T2 + c((cov_MZ_P15),(cov_MZ_P15),(cov_DZ_P15),(cov_DZ_P15),(cov_DZ_P15))*A_P5_T2 + c((cov_MZ_P16),(cov_MZ_P16),(cov_DZ_P16),(cov_DZ_P16),(cov_DZ_P16))*A_P6_T2
    A_P2_T1 ~~ c((cov_MZ_P12),(cov_MZ_P12),(cov_DZ_P12),(cov_DZ_P12),(cov_DZ_P12))*A_P1_T2 + c((cov_MZ_P23),(cov_MZ_P23),(cov_DZ_P23),(cov_DZ_P23),(cov_DZ_P23))*A_P3_T2 + c((cov_MZ_P24),(cov_MZ_P24),(cov_DZ_P24),(cov_DZ_P24),(cov_DZ_P24))*A_P4_T2 + c((cov_MZ_P25),(cov_MZ_P25),(cov_DZ_P25),(cov_DZ_P25),(cov_DZ_P25))*A_P5_T2 + c((cov_MZ_P26),(cov_MZ_P26),(cov_DZ_P26),(cov_DZ_P26),(cov_DZ_P26))*A_P6_T2
    A_P3_T1 ~~ c((cov_MZ_P13),(cov_MZ_P13),(cov_DZ_P13),(cov_DZ_P13),(cov_DZ_P13))*A_P1_T2 + c((cov_MZ_P23),(cov_MZ_P23),(cov_DZ_P23),(cov_DZ_P23),(cov_DZ_P23))*A_P2_T2 + c((cov_MZ_P34),(cov_MZ_P34),(cov_DZ_P34),(cov_DZ_P34),(cov_DZ_P34))*A_P4_T2 + c((cov_MZ_P35),(cov_MZ_P35),(cov_DZ_P35),(cov_DZ_P35),(cov_DZ_P35))*A_P5_T2 + c((cov_MZ_P36),(cov_MZ_P36),(cov_DZ_P36),(cov_DZ_P36),(cov_DZ_P36))*A_P6_T2
    A_P4_T1 ~~ c((cov_MZ_P14),(cov_MZ_P14),(cov_DZ_P14),(cov_DZ_P14),(cov_DZ_P14))*A_P1_T2 + c((cov_MZ_P24),(cov_MZ_P24),(cov_DZ_P24),(cov_DZ_P24),(cov_DZ_P24))*A_P2_T2 + c((cov_MZ_P34),(cov_MZ_P34),(cov_DZ_P34),(cov_DZ_P34),(cov_DZ_P34))*A_P3_T2 + c((cov_MZ_P45),(cov_MZ_P45),(cov_DZ_P45),(cov_DZ_P45),(cov_DZ_P45))*A_P5_T2 + c((cov_MZ_P46),(cov_MZ_P46),(cov_DZ_P46),(cov_DZ_P46),(cov_DZ_P46))*A_P6_T2
    A_P5_T1 ~~ c((cov_MZ_P15),(cov_MZ_P15),(cov_DZ_P15),(cov_DZ_P15),(cov_DZ_P15))*A_P1_T2 + c((cov_MZ_P25),(cov_MZ_P25),(cov_DZ_P25),(cov_DZ_P25),(cov_DZ_P25))*A_P2_T2 + c((cov_MZ_P35),(cov_MZ_P35),(cov_DZ_P35),(cov_DZ_P35),(cov_DZ_P35))*A_P3_T2 + c((cov_MZ_P45),(cov_MZ_P45),(cov_DZ_P45),(cov_DZ_P45),(cov_DZ_P45))*A_P4_T2 + c((cov_MZ_P56),(cov_MZ_P56),(cov_DZ_P56),(cov_DZ_P56),(cov_DZ_P46))*A_P6_T2
    A_P6_T1 ~~ c((cov_MZ_P16),(cov_MZ_P16),(cov_DZ_P16),(cov_DZ_P16),(cov_DZ_P16))*A_P1_T2 + c((cov_MZ_P26),(cov_MZ_P26),(cov_DZ_P26),(cov_DZ_P26),(cov_DZ_P26))*A_P2_T2 + c((cov_MZ_P36),(cov_MZ_P36),(cov_DZ_P36),(cov_DZ_P36),(cov_DZ_P36))*A_P3_T2 + c((cov_MZ_P46),(cov_MZ_P46),(cov_DZ_P46),(cov_DZ_P46),(cov_DZ_P46))*A_P4_T2 + c((cov_MZ_P56),(cov_MZ_P56),(cov_DZ_P56),(cov_DZ_P56),(cov_DZ_P56))*A_P5_T2
   
    
    E_P1_T1 ~~ 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2 + 0*E_P6_T2
    E_P2_T1 ~~ 0*E_P1_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2 + 0*E_P6_T2
    E_P3_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P4_T2 + 0*E_P5_T2 + 0*E_P6_T2
    E_P4_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P5_T2 + 0*E_P6_T2
    E_P5_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P6_T2
    E_P6_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 + 0*E_P3_T2 + 0*E_P4_T2 + 0*E_P5_T2
    
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
    
    A_P6_T1 ~~ 0*E_P6_T1 + 0*E_P6_T2
    A_P6_T2 ~~ 0*E_P6_T1 + 0*E_P6_T2
    

    #cross-trait A*G
    A_P1_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2 
    A_P1_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2 

    A_P2_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2 
    A_P2_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2 
    
    A_P3_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2  
    A_P3_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2  
    
    A_P4_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2  
    A_P4_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P5_T1 + 0*E_P5_T2 + 0*E_P6_T1 + 0*E_P6_T2  
    
    A_P5_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P6_T1 + 0*E_P6_T2  
    A_P5_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P6_T1 + 0*E_P6_T2 
    
    A_P6_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2  
    A_P6_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 + 0*E_P4_T1 + 0*E_P4_T2 + 0*E_P5_T1 + 0*E_P5_T2 
    
    #estimates
    #phenotypic variance
    varP1 := (varA_P1+varE_P1)
    varP2 := (varA_P2+varE_P2)
    varP3 := (varA_P3+varE_P3)
    varP4 := (varA_P4+varE_P4)
    varP5 := (varA_P5+varE_P5)
    varP6 := (varA_P6+varE_P6)
    
    #compute genetic covariance (bivariate heritabilities)
    bh2_P12 := covA_P12/ (covA_P12 + covE_P12)
    bh2_P13 := covA_P13/ (covA_P13 + covE_P13)
    bh2_P14 := covA_P14/ (covA_P14 + covE_P14)
    bh2_P15 := covA_P15/ (covA_P15 + covE_P15)
    bh2_P16 := covA_P15/ (covA_P16 + covE_P16)
    bh2_P23 := covA_P23/ (covA_P23 + covE_P23)
    bh2_P24 := covA_P24/ (covA_P24 + covE_P24)
    bh2_P25 := covA_P25/ (covA_P25 + covE_P25)
    bh2_P26 := covA_P26/ (covA_P26 + covE_P26)
    bh2_P34 := covA_P34/ (covA_P34 + covE_P34)
    bh2_P35 := covA_P35/ (covA_P35 + covE_P35)
    bh2_P36 := covA_P36/ (covA_P36 + covE_P36)
    bh2_P45 := covA_P45/ (covA_P45 + covE_P45)
    bh2_P46 := covA_P46/ (covA_P45 + covE_P46)
    bh2_P56 := covA_P56/ (covA_P56 + covE_P56)
    
    #compute genetic correlations
    sdA_P1 := sqrt(varA_P1)
    sdA_P2 := sqrt(varA_P2)
    sdA_P3 := sqrt(varA_P3)
    sdA_P4 := sqrt(varA_P4)
    sdA_P5 := sqrt(varA_P5)
    sdA_P6 := sqrt(varA_P6)
    
    rA_P12 := covA_P12/(sdA_P1*sdA_P2)
    rA_P13 := covA_P13/(sdA_P1*sdA_P3)
    rA_P14 := covA_P14/(sdA_P1*sdA_P4)
    rA_P15 := covA_P15/(sdA_P1*sdA_P5)
    rA_P16 := covA_P16/(sdA_P1*sdA_P6)
    rA_P23 := covA_P23/(sdA_P2*sdA_P3)
    rA_P24 := covA_P24/(sdA_P2*sdA_P4)
    rA_P25 := covA_P25/(sdA_P2*sdA_P5)
    rA_P26 := covA_P26/(sdA_P2*sdA_P6)
    rA_P34 := covA_P34/(sdA_P3*sdA_P4)
    rA_P35 := covA_P35/(sdA_P3*sdA_P5)
    rA_P36 := covA_P36/(sdA_P3*sdA_P6)
    rA_P45 := covA_P45/(sdA_P4*sdA_P4)
    rA_P46 := covA_P46/(sdA_P4*sdA_P6)
    rA_P56 := covA_P56/(sdA_P5*sdA_P6)
    
    #compute environmental correlations
    sdE_P1 := sqrt(varE_P1)
    sdE_P2 := sqrt(varE_P2)
    sdE_P3 := sqrt(varE_P3)
    sdE_P4 := sqrt(varE_P4)
    sdE_P5 := sqrt(varE_P5)
    sdE_P6 := sqrt(varE_P6)
    
    rE_P12 := covE_P12/(sdE_P1*sdE_P2)
    rE_P13 := covE_P13/(sdE_P1*sdE_P3)
    rE_P14 := covE_P14/(sdE_P1*sdE_P4)
    rE_P15 := covE_P15/(sdE_P1*sdE_P5)
    rE_P16 := covE_P16/(sdE_P1*sdE_P6)
    rE_P23 := covE_P23/(sdE_P2*sdE_P3)
    rE_P24 := covE_P24/(sdE_P2*sdE_P4)
    rE_P25 := covE_P25/(sdE_P2*sdE_P5)
    rE_P26 := covE_P26/(sdE_P2*sdE_P6)
    rE_P34 := covE_P34/(sdE_P3*sdE_P4)
    rE_P35 := covE_P35/(sdE_P3*sdE_P5)
    rE_P36 := covE_P36/(sdE_P3*sdE_P6)
    rE_P45 := covE_P45/(sdE_P4*sdE_P5)
    rE_P46 := covE_P46/(sdE_P4*sdE_P6)
    rE_P56 := covE_P56/(sdE_P5*sdE_P6)
    
    #twin-h2
    h2_P1 := varA_P1/ varP1
    h2_P2 := varA_P2/ varP2
    h2_P3 := varA_P3/ varP3
    h2_P4 := varA_P4/ varP4
    h2_P5 := varA_P5/ varP5
    h2_P6 := varA_P6/ varP6
    
    #covA difference
    dcovA_P123 := covA_P12 - covA_P13
    dcovA_P124 := covA_P12 - covA_P14
    dcovA_P125 := covA_P12 - covA_P15
    dcovA_P126 := covA_P12 - covA_P16
    dcovA_P134 := covA_P13 - covA_P14
    dcovA_P135 := covA_P13 - covA_P15
    dcovA_P136 := covA_P13 - covA_P16
    dcovA_P145 := covA_P14 - covA_P15
    dcovA_P146 := covA_P14 - covA_P16
    dcovA_P156 := covA_P15 - covA_P16
    
    #transform to rA zA
    zA_P12 := 0.5 * log((1 + rA_P12)/(1 - rA_P12))
    zA_P13 := 0.5 * log((1 + rA_P13)/(1 - rA_P13))
    zA_P14 := 0.5 * log((1 + rA_P14)/(1 - rA_P14))
    zA_P15 := 0.5 * log((1 + rA_P15)/(1 - rA_P15))
    zA_P16 := 0.5 * log((1 + rA_P16)/(1 - rA_P16))

    #zA differences (note that correlations are dependent)
    dzA_P123 := zA_P12 - zA_P13
    dzA_P124 := zA_P12 - zA_P14
    dzA_P125 := zA_P12 - zA_P15
    dzA_P126 := zA_P12 - zA_P16
    dzA_P134 := zA_P13 - zA_P14
    dzA_P135 := zA_P13 - zA_P15
    dzA_P136 := zA_P13 - zA_P16
    dzA_P145 := zA_P14 - zA_P15
    dzA_P146 := zA_P14 - zA_P16
    dzA_P156 := zA_P15 - zA_P16
    
    #covE difference
    dcovE_P123 := covE_P12 - covE_P13
    dcovE_P124 := covE_P12 - covE_P14
    dcovE_P125 := covE_P12 - covE_P15
    dcovE_P126 := covE_P12 - covE_P16
    dcovE_P134 := covE_P13 - covE_P14
    dcovE_P135 := covE_P13 - covE_P15
    dcovE_P136 := covE_P13 - covE_P16
    dcovE_P145 := covE_P14 - covE_P15
    dcovE_P146 := covE_P14 - covE_P16
    dcovE_P156 := covE_P15 - covE_P16
    
    #transform to rE zE
    zE_P12 := 0.5 * log((1 + rE_P12)/(1 - rE_P12))
    zE_P13 := 0.5 * log((1 + rE_P13)/(1 - rE_P13))
    zE_P14 := 0.5 * log((1 + rE_P14)/(1 - rE_P14))
    zE_P15 := 0.5 * log((1 + rE_P15)/(1 - rE_P15))
    zE_P16 := 0.5 * log((1 + rE_P16)/(1 - rE_P16))

    #zE differences (note that correlations are dependent)
    dzE_P123 := zE_P12 - zE_P13
    dzE_P124 := zE_P12 - zE_P14
    dzE_P125 := zE_P12 - zE_P15
    dzE_P126 := zE_P12 - zE_P16
    dzE_P134 := zE_P13 - zE_P14
    dzE_P135 := zE_P13 - zE_P15
    dzE_P136 := zE_P13 - zE_P16
    dzE_P145 := zE_P14 - zE_P15
    dzE_P146 := zE_P14 - zE_P16
    dzE_P156 := zE_P15 - zE_P16
  "

#  Test significance of difference between genetic correlations by constraining A covariances to be equal
constrain_AP123 = " covA_P12 == covA_P13 "
constrain_AP124 = " covA_P12 == covA_P14 "
constrain_AP125 = " covA_P12 == covA_P15 "
constrain_AP126 = " covA_P12 == covA_P16 "
constrain_AP134 = " covA_P13 == covA_P14 "
constrain_AP135 = " covA_P13 == covA_P15 "
constrain_AP136 = " covA_P13 == covA_P16 "
constrain_AP145 = " covA_P14 == covA_P15 "
constrain_AP146 = " covA_P14 == covA_P16 "
constrain_AP156 = " covA_P15 == covA_P16 "

#  Test significance of difference between genetic correlations by constraining E covariances to be equal
constrain_EP123 = " covE_P12 == covE_P13 "
constrain_EP124 = " covE_P12 == covE_P14 "
constrain_EP125 = " covE_P12 == covE_P15 "
constrain_EP126 = " covE_P12 == covE_P16 "
constrain_EP134 = " covE_P13 == covE_P14 "
constrain_EP135 = " covE_P13 == covE_P15 "
constrain_EP136 = " covE_P13 == covE_P16 "
constrain_EP145 = " covE_P14 == covE_P15 "
constrain_EP146 = " covE_P14 == covE_P16 "
constrain_EP156 = " covE_P15 == covE_P16 "

#  FOR NOW THIS MODEL IS ONLY SPECIFIED FOR THE CFS