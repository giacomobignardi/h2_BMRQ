#Description: Install packages in the VDI
#this script contains all the R packages needed to run the main analyses
#Program: R_packages ------------------------------------------------------------------------------------------------------------------------------
#set path where to store the packages
.libPaths('W:\\C4_Behavior_Genectics_Unit\\Giacomo Bignardi\\2023_BMRQ\\04_packages')
#get packages from KI repo
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org/')
install.packages("lavaan", dependencies = TRUE)
install.packages("tidyverse")
install.packages("rstatix")
install.packages("psych")
install.packages("patchwork")
install.packages("openxlsx")



