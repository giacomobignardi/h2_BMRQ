#ADE MODEL####
# structural equation model specification for 5 groups and one covariate:
# 1:monozygotic female
# 2:monozygotic male
# 3:dizygotic female
# 4:dizygotic male
# 5:dizygotic opposite-sex
# based on SAT results gropus are constrained to 2 (5 gropus constrained)

#ADE MODEL####
#full ADE model with covariate
ADE_model <-"
 #estimate intercepts (means) for males and females separately adjusting for age as covariate
P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
#constrain covariate variances to be equal across groups
age ~~ var_age*age
#ADE components
A_T1=~ P_T1
A_T2=~ P_T2
D_T1 =~ P_T1
D_T2 =~ P_T2
E_T1 =~ P_T1
E_T2 =~ P_T2
#ADE variances
A_T1 ~~ varA*A_T1
A_T2 ~~ varA*A_T2
D_T1 ~~ varD*D_T1
D_T2 ~~ varD*D_T2
E_T1 ~~ varE*E_T1
E_T2 ~~ varE*E_T2
#fix remaining P variance to 0
P_T1~~0*P_T1
P_T2~~0*P_T2
#CTD twin informed constraints
covMZA == varA 
covDZA == .5*varA 
covMZD == varD
covDZD == .25*varD
#covariances
#covariances cross-twin wihtin-trait
    A_T1 ~~ c((covMZA),(covMZA),(covDZA),(covDZA),(covDZA))*A_T2
    D_T1 ~~ c((covMZD),(covMZD),(covDZD),(covDZD),(covDZD))*D_T2
    E_T1 ~~ 0*E_T2
#set every other covariance to zero
A_T1 ~~ 0*E_T1 + 0*E_T2 + 0*D_T1 + 0*D_T2
A_T2 ~~ 0*E_T1 + 0*E_T2 + 0*D_T1 + 0*D_T2
D_T1 ~~ 0*E_T1 + 0*E_T2
D_T2 ~~ 0*E_T1 + 0*E_T2

#calculate summary output
varP := (varA+varD+varE)
sdA := sqrt(varA)
sdD := sqrt(varD)
sdE := sqrt(varE)

#narrow-sense twin-h2
h2_P := varA / varP
A := varA / varP
D := varD / varP
E := varE / varP
"
#AE MODEL####
#constained AE model with covariate
AE_model <-"
#estimate intercepts (means) for males and females separately adjusting for age as covariate
P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
#constrain covariate variances to be equal across groups
age ~~ var_age*age
#ADE components
A_T1=~ P_T1
A_T2=~ P_T2
E_T1 =~ P_T1
E_T2 =~ P_T2
#ADE variances
A_T1 ~~ varA*A_T1
A_T2 ~~ varA*A_T2
E_T1 ~~ varE*E_T1
E_T2 ~~ varE*E_T2
#fix remaining P variance to 0
P_T1~~0*P_T1
P_T2~~0*P_T2
#CTD twin informed constraints
covMZA == varA 
covDZA == .5*varA 
#covariances
#covariances cross-twin wihtin-trait
    A_T1 ~~ c((covMZA),(covMZA),(covDZA),(covDZA),(covDZA))*A_T2
    E_T1 ~~ 0*E_T2
#set every other covariance to zero
A_T1 ~~ 0*E_T1 + 0*E_T2 
A_T2 ~~ 0*E_T1 + 0*E_T2

#calculate summary output
varP := (varA+varE)
sdA := sqrt(varA)
sdE := sqrt(varE)

#narrow-sense twin-h2
h2_P := varA / varP
A := varA / varP
E := varE / varP
"
#E MODEL####
#constained E model with covariate
E_model <-"
#estimate intercepts (means) for males and females separately adjusting for age as covariate
P_T1 ~ c(mean_f,mean_m,mean_f,mean_m,mean_f)*1 + b1*age
P_T2 ~ c(mean_f,mean_m,mean_f,mean_m,mean_m)*1 + b1*age
#constrain covariate variances to be equal across groups
age ~~ var_age*age
#ADE components
E_T1 =~ P_T1
E_T2 =~ P_T2
#ADE variances
E_T1 ~~ varE*E_T1
E_T2 ~~ varE*E_T2
#fix remaining P variance to 0
P_T1~~0*P_T1
P_T2~~0*P_T2
#covariances
#covariances cross-twin wihtin-trait
E_T1 ~~ 0*E_T2
#calculate summary output
varP := varE
sdE := sqrt(varE)
E := varE / varP
"