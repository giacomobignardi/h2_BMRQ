#Trivariate Sequential####
# structural equation model specification for 5 groups (constrained to 2) and one covariate:
# 1:monozygotic female
# 2:monozygotic male
# 3:dizygotic female
# 4:dizygotic male
# 5:dizygotic opposite-sex

Seq_AE_mod <-"
    #create measurement model for twin one and twin two
    #mean (allow means to be different across sex)
    P1_T1~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_f_P1)*1 + b1_P1*age
    P1_T2~c(mean_f_P1, mean_m_P1, mean_f_P1, mean_m_P1, mean_m_P1)*1 + b1_P1*age
    P2_T1~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_f_P2)*1 + b1_P2*age
    P2_T2~c(mean_f_P2, mean_m_P2, mean_f_P2, mean_m_P2, mean_m_P2)*1 + b1_P2*age
    P3_T1~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_f_P3)*1 + b1_P3*age
    P3_T2~c(mean_f_P3, mean_m_P3, mean_f_P3, mean_m_P3, mean_m_P3)*1 + b1_P3*age
    
    #constrain covariate variances to be equal across groups
    age ~~ var_age*age
    

    #AE components
    A_P1_T1=~ P1_T1
    A_P1_T2=~ P1_T2
    A_P2_T1=~ P2_T1
    A_P2_T2=~ P2_T2
    A_P3_T1=~ P3_T1
    A_P3_T2=~ P3_T2

    E_P1_T1 =~ P1_T1
    E_P1_T2 =~ P1_T2
    E_P2_T1 =~ P2_T1
    E_P2_T2 =~ P2_T2
    E_P3_T1 =~ P3_T1
    E_P3_T2 =~ P3_T2

    # variances
    A_P1_T1 ~~ varA_P1*A_P1_T1
    A_P1_T2 ~~ varA_P1*A_P1_T2
    A_P2_T1 ~~ varA_P2*A_P2_T1
    A_P2_T2 ~~ varA_P2*A_P2_T2
    A_P3_T1 ~~ varA_P3*A_P3_T1
    A_P3_T2 ~~ varA_P3*A_P3_T2

    E_P1_T1~~varE_P1*E_P1_T1
    E_P1_T2~~varE_P1*E_P1_T2
    E_P2_T1~~varE_P2*E_P2_T1
    E_P2_T2~~varE_P2*E_P2_T2
    E_P3_T1~~varE_P3*E_P3_T1
    E_P3_T2~~varE_P3*E_P3_T2

    P1_T1~~0*P1_T1
    P1_T2~~0*P1_T2
    P2_T1~~0*P2_T1
    P2_T2~~0*P2_T2
    P3_T1~~0*P3_T1
    P3_T2~~0*P3_T2

    #constraints (twin pair covariances, p = phenotype)
    covMZ_p1 == varA_P1
    covMZ_p2 == varA_P2
    covMZ_p3 == varA_P3

    covDZ_p1 == .5*varA_P1
    covDZ_p2 == .5*varA_P2
    covDZ_p3 == .5*varA_P3

    # covariances
    A_P1_T1 ~~ c((covMZ_p1),(covMZ_p1),(covDZ_p1),(covDZ_p1),(covDZ_p1))*A_P1_T2
    A_P2_T1 ~~ c((covMZ_p2),(covMZ_p2),(covDZ_p2),(covDZ_p2),(covDZ_p2))*A_P2_T2
    A_P3_T1 ~~ c((covMZ_p3),(covMZ_p3),(covDZ_p3),(covDZ_p3),(covDZ_p3))*A_P3_T2

    #set every other covariance to 0
    #within trait E
    E_P1_T1 ~~ 0*E_P1_T2
    E_P2_T1 ~~ 0*E_P2_T2
    E_P3_T1 ~~ 0*E_P3_T2

    #within trait A*G
    A_P1_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2
    A_P1_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2

    A_P2_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2
    A_P2_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2

    A_P3_T1 ~~ 0*E_P3_T1 + 0*E_P3_T2
    A_P3_T2 ~~ 0*E_P3_T1 + 0*E_P3_T2

    #cross trait A*G
    A_P1_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2
    A_P1_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2

    A_P2_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2
    A_P2_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2

    A_P3_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2
    A_P3_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2

    #cross trait A
    A_P1_T1 ~~ 0*A_P2_T1 + 0*A_P2_T2 + 0*A_P3_T1 + 0*A_P3_T2
    A_P1_T2 ~~ 0*A_P2_T1 + 0*A_P2_T2 + 0*A_P3_T1 + 0*A_P3_T2
    A_P3_T1 ~~ 0*A_P2_T1 + 0*A_P2_T2
    A_P3_T2 ~~ 0*A_P2_T1 + 0*A_P2_T2

    #cross trait E
    E_P1_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2
    E_P1_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2
    E_P3_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2
    E_P3_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2

    #path
    P2_T1 ~ a12*A_P1_T1
    P2_T2 ~ a12*A_P1_T2
    P2_T1 ~ e12*E_P1_T1
    P2_T2 ~ e12*E_P1_T2
    P3_T1 ~ a13*A_P1_T1 + a23*A_P2_T1
    P3_T2 ~ a13*A_P1_T2 + a23*A_P2_T2
    P3_T1 ~ e13*E_P1_T1 + e23*E_P2_T1
    P3_T2 ~ e13*E_P1_T2 + e23*E_P2_T2

    # meaningful if not Cholesky (standardisedsolution)
    #variance 
    varP1 := varA_P1 + varE_P1
    varP2 := varA_P2 + varE_P2 + (varA_P1*a12^2) + (varE_P1*e12^2)
    varP3 := varA_P3 + varE_P3 + (varA_P1*a13^2) + (varE_P1*e13^2) + (varA_P2*a23^2) + (varE_P2*e23^2)
    
    #variance components
    h2_P1 := varA_P1/varP1
    h2_P2 := varA_P2/varP2
    h2_P3 := varA_P3/varP3
    e1_P3 := varE_P1/varP3
    e2_P3 := varE_P2/varP3
    e3_P3 := varE_P3/varP3
    
    #proportion of h2
    varA31ut := (varA_P1*a13^2)/(varA_P3 + (varA_P1*a13^2) + (varA_P2*a23^2))
    varA32ut := (varA_P2*a23^2)/(varA_P3 + (varA_P1*a13^2) + (varA_P2*a23^2))
    varA33ut := varA_P3/(varA_P3 + (varA_P1*a13^2) + (varA_P2*a23^2))

    #R2
    R2a12 := (varA_P1*a12^2)/varP2
    R2a13 := (varA_P1*a13^2)/varP3
    R2a23 := (varA_P2*a23^2)/varP3
    R2e12 := (varE_P1*e12^2)/varP2
    R2e13 := (varE_P1*e13^2)/varP3
    R2e23 := (varE_P2*e23^2)/varP3
    R2ae123:= ((varA_P1*a13^2) + (varE_P1*e13^2) + (varA_P2*a23^2) + (varE_P2*e23^2))/varP3
  "

#test path significance by LRT
#path from A1
lrt_a12 <-"a12==0"
lrt_a13 <-"a13==0"
#path from A2
lrt_a23 <-"a23==0"
#path from E1
lrt_e12 <-"e12==0"
lrt_e13 <-"e13==0"
#path from E2
lrt_e23 <-"e23==0"