#Description: Install packages in the VDI
#this script contains all the extra R packages needed to run the analyses (including supplementary analyses)
#Program: R_packages ------------------------------------------------------------------------------------------------------------------------------
#set path where to store the packages
.libPaths('W:\\XXX\04_packages')
#get packages from KI repo
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org/')
install.packages("lavaan", dependencies = TRUE)
install.packages("tidyverse")
install.packages("rstatix")
install.packages("psych")
install.packages("patchwork")
install.packages("wesanderson")
install.packages("OpenMx")
install.packages("ggnewscale")
install.packages('umx')



