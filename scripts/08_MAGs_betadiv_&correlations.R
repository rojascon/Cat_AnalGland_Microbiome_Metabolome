################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 22 Apr 2023
#                       Last updated: 22 April 2023
#
################################################################################

## CODE FOR: A) determining whether MAG relative abundances vary with host
#               characteristics (age, diet, obesity, environment, periodontal disease)

##           B) determining whether MAG abundance profiles are correlated with
#                   metabolome profiles extracted during solid-phase or
#                     after liquid derivatization using GC-MS


################################################################################
#             1.  Load MAG abundances and sample metadata
################################################################################\

# MAG abundances
coverm=read.table("data/00_MAG_Abun_SPADES.txt", sep="\t", header=T,row.names=1);

# sample metadata
load("data/01_metadata_formatted.Rdata");

# make a vector of your sample names 
meta2=meta[meta$swab_type=="Anal.Sac",];
coverm=coverm[,meta2$sampleID];

# see which MAGs are the most abundant
mydf=coverm[ order(rowMeans(coverm), decreasing = T),]
rowMeans(mydf)[1:20];


################################################################################
#           2. Generate Bray-Curtis distances based on 
#           MAG abundances across samples and run PERMANOVA
################################################################################

# relative abundances for Bray-Curtis
kspe=coverm;
bray<-apply(kspe, 2, function(i) (i/sum(i))*100);

# presence absence for Jaccard
jac=(kspe>0)*1;

# clr (center log ratio) for Aitchison distance
clra=data.frame(clr(kspe));

#generate distance matrices
ps1<- phyloseq(otu_table(bray, taxa_are_rows=TRUE));
bray.dist=phyloseq::distance(ps1, method="bray");

ps2<- phyloseq(otu_table(jac, taxa_are_rows=TRUE));
jac.dist=phyloseq::distance(ps2, method="jaccard",binary=T);

ps3<- phyloseq(otu_table(clra, taxa_are_rows=TRUE));
clr.dist=phyloseq::distance(ps3, method="euclidean");

# run PERMANOVAs
mydist=list(jac.dist, bray.dist,clr.dist);
names=c("Jaccard","Bray-Curtis","Aitchison");
met=c("jaccard","bray","aitchison");

for(i in 1:3)
{
  print(paste("PERMANOVA test, MAGs, N=23, using:", names[i]));
  print(adonis2(mydist[[i]]~
                  age_yrs+
                  obese+
                  environment+
                  diet_bi+
                  Periodontal,
                data=meta2,
                method = met[i],
                by="margin",
                permutations = 999));
};


################################################################################
#           3. Correlate MAG abundances profiles with metabolome profiles
################################################################################

# using jaccard and euclidean distances for volatile data 
# that were calculated in script 04_metabolome_microbiome_correlations.R

# Mantel tests (volatile data solid phase ~ microbiome data);
load("data/04_metabolite_solid_distances.Rdata");
mantel(vjac.dist,jac.dist, method="spearman", permutations=999);
mantel(veuc.dist,bray.dist, method="spearman", permutations=999);
mantel(veuc.dist,clr.dist, method="spearman", permutations=999);

# Mantel tests (volatile data liquid deriv. ~ microbiome data);
load("data/04_metabolite_liquid_distances.Rdata");
mantel(vjac.dist,jac.dist, method="spearman", permutations=999);
mantel(veuc.dist,bray.dist, method="spearman", permutations=999);
mantel(veuc.dist,clr.dist, method="spearman", permutations=999);

