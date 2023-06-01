################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 22 April 2023
#                       Last updated: 23 April 2023
#
################################################################################

## CODE FOR:  determining whether the microbiome in the anal gland is distinct
#             from the microbiome in the perianal region in terms of their:
#             a) microbiome compositions (Kraken Genus-abundances)
#             b) microbiome functions (Anvio COG & KEGG abundances)
#             c) MAG abundances (metagenome-assembled genomes);


source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1.  Load sample metadata, Kraken microbiome data, 
#                   COG & KEGG abundances, and MAG abundances
################################################################################

# sample metadata
load("data/01_metadata_formatted.Rdata");
meta2=meta[meta$name %in% c("Cannoli","Cheddar","Cooper",
                            "Gabby","Goose","Virginia"),];

# Kraken abundances
load("data/01_bracken_Genus_nocontam.Rdata");

# MAG abundances
coverm=read.table("data/00_MAG_Abun_SPADES.txt", 
                  sep="\t", header=T,row.names=1);

# pick one of the following: 
      # COG and KEGG (broad pathways)
load("data/06_COGpath_spades.Rdata"); #70 rows
load("data/06_KOModule_abundances_spades.Rdata"); #192 rows

      # COG and KEGG (specific genes)
load("data/06_COGfun_spades.Rdata"); #4150 rows
load("data/06_KOfam_abundances_spades.Rdata"); #5509 rows


################################################################################
#           2. Kraken microbiome abundances ~ body site
################################################################################

# aggregate microbiome taxonomic data by Genus
kspe=krak[,which(names(krak) 
                 %in% c(meta2$sampleID, "Genus"))];

colnames(kspe)[1]="taxa";
kspe=aggregate(.~taxa,kspe, sum);
rownames(kspe)=kspe$taxa; kspe$taxa=NULL;
kspe=kspe[rowSums(kspe)>0,];
kspe=kspe[,order(colnames(kspe))];

# relative abundances for Bray-Curtis
bray<-apply(kspe, 2, function(i) (i/sum(i))*100);

# presence absence for Jaccard
jac=(kspe>0)*1;

# clr (center log ratio) for Aitchison distance ("compositions")
clra=data.frame(compositions::clr(kspe));

# generate distance matrices
ps1<- phyloseq(otu_table(bray, taxa_are_rows=TRUE));
bray.dist=phyloseq::distance(ps1, method="bray");

ps2<- phyloseq(otu_table(jac, taxa_are_rows=TRUE));
jac.dist=phyloseq::distance(ps2, method="jaccard",binary=T);

ps3<- phyloseq(otu_table(clra, taxa_are_rows=TRUE));
clr.dist=phyloseq::distance(ps3, method="euclidean");

# run PERMANOVAs using for loop
mydist=list(jac.dist, bray.dist,clr.dist);
names=c("Jaccard","Bray-Curtis","Aitchison");
met=c("jaccard","bray","aitchison");

for(i in 1:3)
{
  print(paste("PERMANOVA test, perianal vs. microbiome, N=12, using:", names[i]));
  print(adonis2(mydist[[i]]~
                  swab_type+
                  name,
                data=meta2,
                method = met[i],
                by="margin",
                permutations = 999));
};


################################################################################
#           2. Anvio microbiome functions ~ body site
################################################################################

# clean up COG and KEGG tables
cogabun=cogabun[,order(colnames(cogabun))];
keggabun=keggabun[,order(colnames(keggabun))];

rownames(cogabun)=cogabun$gf; cogabun$gf=NULL;
rownames(keggabun)=keggabun$gf; keggabun$gf=NULL;

cogabun=cogabun[rowSums(cogabun)>0,];
keggabun=keggabun[rowSums(keggabun)>0,];

# retain only samples from the six cats
cogabun=cogabun[,meta2$sampleID];
keggabun=keggabun[,meta2$sampleID];

# calculate proportions for Bray distance
cbray<-apply(cogabun, 2, function(i) (i/sum(i))*100);
kbray<-apply(keggabun, 2, function(i) (i/sum(i))*100);

# convert to presence absence for Jaccard distance
cjac=(cogabun>0)*1;
kjac=(keggabun>0)*1;

# clr (center log ratio) tranformation for Aitchison distance 
cclra=data.frame(clr(cogabun));
kclra=data.frame(clr(keggabun));

# make microbiome functions distance matrices
ps1<- phyloseq(otu_table(cjac, taxa_are_rows=TRUE));
cjac.dist=phyloseq::distance(ps1, method="jaccard",binary=T);
ps2<- phyloseq(otu_table(kjac, taxa_are_rows=TRUE));
kjac.dist=phyloseq::distance(ps2, method="jaccard",binary=T);

ps3<- phyloseq(otu_table(cbray, taxa_are_rows=TRUE));
cbray.dist=phyloseq::distance(ps3, method="bray");
ps4<- phyloseq(otu_table(kbray, taxa_are_rows=TRUE));
kbray.dist=phyloseq::distance(ps4, method="bray");

ps5<- phyloseq(otu_table(cclra, taxa_are_rows=TRUE));
cclr.dist=phyloseq::distance(ps5, method="euclidean");
ps6<- phyloseq(otu_table(kclra, taxa_are_rows=TRUE));
kclr.dist=phyloseq::distance(ps6, method="euclidean");

# run PERMANOVAs using for loop
mydist=list(cjac.dist,kjac.dist, cbray.dist, 
            kbray.dist, cclr.dist, kclr.dist);
names=c("CJaccard","KJaccard","CBray-Curtis","KBray-Curtis","CAitchison","KAitchison");
met=c("jaccard","jaccard","bray","bray","aitchison","aitchison");

for(i in 1:6)
{
  print(paste("PERMANOVA test, microbiome functions, N=12, using:", names[i]));
  print(adonis2(mydist[[i]]~
                  swab_type+
                  name,
                data=meta2,
                method = met[i],
                by="margin",
                permutations = 999));
};


################################################################################
#           3. MAG abundances ~ body site
################################################################################

# retain samples from the six cats 
coverm=coverm[,meta2$sampleID];

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

# run PERMANOVAs using for loop
mydist=list(jac.dist, bray.dist,clr.dist);
names=c("Jaccard","Bray-Curtis","Aitchison");
met=c("jaccard","bray","aitchison");

for(i in 1:3)
{
  print(paste("PERMANOVA test, MAGs perianal vs. microbiome, N=12, using:", names[i]));
  print(adonis2(mydist[[i]]~
                  swab_type+
                  name,
                data=meta2,
                method = met[i],
                by="margin",
                permutations = 999));
};