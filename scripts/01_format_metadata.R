################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 17 July 2022
#                       Last updated: 21 April 2023
#
################################################################################

## CODE FOR: formatting sample metadata factors

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1.  Load sample metadata and only retain samples of interest                
################################################################################

#load metadata
meta=read.csv("data/00_sample_metadata.csv", stringsAsFactors = F);

################################################################################
#             2.  Format metadata factors                
################################################################################

meta$DOB=as.Date(meta$DOB, format="%d-%b-%y");
meta$sample_date=as.Date(meta$sample_date, format="%d-%b-%y");
meta$age_cat=factor(meta$age_cat, levels=c("2-5","6-8","9",">10"));
meta$obese=factor(meta$obese);
meta$environment=factor(meta$environment);
meta$diet_bi=factor(meta$diet_bi, levels=c("Dry","Other"));
meta$Periodontal=factor(meta$Periodontal, levels=c("Y","N"));

meta=meta[order(meta$sampleID),];

################################################################################
#             3. save formatted metadata file
################################################################################

save(meta, file="data/01_metadata_formatted.Rdata");

