################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 30 Jan 2023
#                       Last updated: 3 Feb 2023
#
################################################################################

## CODE FOR: removing bacterial species that might be contaminants
#            from the dataset using information obtained from 2 sterile swab samples

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1.  Load taxa abundance table output by Kraken 
#                           and sample metadata
################################################################################

# Species-level abundances
krak=read.csv("data/00_bracken_Genus.csv", header=T);
krak$Genus=NULL; 

# sample metadata
load("data/01_metadata_formatted.Rdata");
rownames(meta)=meta$sampleID; meta$sampleID=NULL;

# convert raw counts to relative abundances
colnames(krak)[1]="taxa";
krak=aggregate(.~taxa,krak, sum);
krak=krak[rowSums(krak[-1])>0,];
rownames(krak)=krak$taxa; krak$taxa=NULL;

# remove perianal samples temporarily
meta2=meta[meta$swab_type!="Perianal",];
krak2=krak[,colnames(krak)%in% rownames(meta2)];


################################################################################
#             2.  Identify potential contaminant bacterial species using
#                                   R decontam 
################################################################################

# combine the 2 tables into a phyloseq object
ps=phyloseq(otu_table(krak, 
                      taxa_are_rows =T),
            sample_data(meta));

# identify potential contaminants using R decontam package
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# using the "prevalence" method & a 0.5 decontam score threshold

sample_data(ps)$is.neg=sample_data(ps)$swab_type == "Blank";
scores <- isContaminant(ps, method="prevalence", neg="is.neg",
                        threshold = 0.5, detailed = TRUE);

# view # of contaminants identified by decontam
table(scores$contaminant);


################################################################################
#             3.  Remove these contaminant bacterial species 
################################################################################

# filter out contaminants with decontam scores <0.5 
contam=scores[scores$contaminant=="TRUE",];
krak=krak[!rownames(krak) %in% rownames(contam),];
krak=krak[rowSums(krak)>0,];

# bring back the 'Genus' column
glab=as.data.frame(str_split_fixed(rownames(krak)," ",2));
glab$Species=rownames(krak); glab$V2=NULL; colnames(glab)[1]="Genus";
krak=cbind(glab,krak);
rownames(krak)=NULL;

# save file
save(krak,file="data/01_bracken_Genus_nocontam.Rdata");
