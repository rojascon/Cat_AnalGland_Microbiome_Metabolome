################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 16 July 2022
#                       Last updated: 31 May 2023
#
################################################################################

## CODE FOR: configuring R workspace and printing R version and package versions
#for reader


################################################################################
#             1.  Configure the workspace for subsequent R project scripts                 
################################################################################

#set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE) ;

#load necessary packages
library(pacman);
pacman::p_load("phyloseq","reshape2","ggplot2","tidyr", "vegan","decontam", 
               "compositions", "gridExtra", "forcats",
               "microbiome");


################################################################################
#             2. Communicate the R version and package versions to reader                 
################################################################################

print("This code was developed with R version 4.3.0");

print("The packages used and their versions were: phyloseq_1.44.0| reshape2_1.4.4|
      ggplot2_3.4.2| tidyr_1.3.0| vegan_2.6-4| decontam_1.20.0| compositions_2.0-6|
      gridExtra_2.3| forcats_1.0.0| microbiome_1.22.0");

