################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 16 Aug 2022
#                       Last updated: 22 April 2023
#
################################################################################

## CODE FOR: A) testing whether anal gland microbiome functional profiles vary with
#            host characteristics (age, diet, obesity, environment, periodontal disease)
#            Cluster of Orthologous Genes (COG); Kyoto Encyclopedia Genes Genomes (KEGG)

#            B) correlating microbiome functional profiles with metabolome profiles
#                 VOCs -- volatile organic compounds
#            C) generating stacked barplots of KEGG and COG pathway abundances

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1. Load KEGG and COG data, and sample metadata               
################################################################################

# pick one of the following either broad or specific: 
    # COG and KEGG (broad pathways)
load("data/06_COGpath_spades.Rdata"); #70 rows
load("data/06_KOModule_abundances_spades.Rdata"); #192 rows

    # COG and KEGG (specific genes)
load("data/06_COGfun_spades.Rdata"); #4150 rows
load("data/06_KOfam_abundances_spades.Rdata"); #5509 rows

# sample metadata
load("data/01_metadata_formatted.Rdata");


################################################################################
#             2. Calculate KEGG and COG relative abundances and do the necessary
#               data transformations or scaling in preparation for calculating
#                               beta-diversity distances
################################################################################

# clean up COG and KEGG tables
cogabun=cogabun[,order(colnames(cogabun))];
keggabun=keggabun[,order(colnames(keggabun))];

rownames(cogabun)=cogabun$gf; cogabun$gf=NULL;
rownames(keggabun)=keggabun$gf; keggabun$gf=NULL;

cogabun=cogabun[rowSums(cogabun)>0,];
keggabun=keggabun[rowSums(keggabun)>0,];

# retain only samples from the anal gland
meta2=meta[meta$swab_type=="Anal.Sac",];
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


################################################################################
#             3. Run microbiome functions PERMANOVAS
#       microbiome functional beta-diversity ~ host factors
################################################################################
      
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
  print(paste("PERMANOVA test, microbiome functions, N=23, using:", names[i]));
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

# make sure to run PERMANOVAs for both specific COG/KEGG genes and broad COG/KEGG pathways

# calculate means of the most abundant gene functions 
mydf=cbray;   # cbray or kbray
mydf=mydf[ order(rowMeans(mydf), decreasing = T),]
rowMeans(mydf)[1:20];


################################################################################
#             3. correlate microbiome functional profiles with metabolome profiles            
################################################################################

# load jaccard and euclidean distances for volatile data 
# calculated in script 04_metabolome_microbiome_correlations.R
# pick one: either solid-phase or liquid derivatization
load("data/04_metabolite_solid_distances.Rdata");
load("data/04_metabolite_liquid_distances.Rdata");

# jaccard VOCs ~ jaccard COG or jaccard KEGG
mantel(vjac.dist,cjac.dist, method="spearman", permutations=999);
mantel(vjac.dist,kjac.dist, method="spearman", permutations=999);

# euclidean VOCs ~ bray COG or bray KEGG
mantel(veuc.dist,cbray.dist, method="spearman", permutations=999);
mantel(veuc.dist,kbray.dist, method="spearman", permutations=999);

# euclidean VOCs ~ euclidean COG or euclidean KEGG
mantel(veuc.dist,cclr.dist, method="spearman", permutations=999);
mantel(veuc.dist,kclr.dist, method="spearman", permutations=999);

# remember to rerun the above 6 Mantel tests for both solid-phase and liquid-phase VOC data
# and for both specific functional genes and broad functional pathways
# so a total of 4 times


################################################################################
#             5. Barplot of KEGG and COG pathway abundances
#                   of pathways involved in odor-compound synthesis
################################################################################

# make dataframe of the pathways and their mean relative abundances
df=data.frame("Predictors"=c("fatty-acid biosynthesis","lipid A biosynthesis", 
                             "aromatic amino acid biosynthesis",
                             "ketone body biosynthesis","cholesterol biosynthesis",
                             "fatty acid biosynthesis","lipid A biosynthesis",
                             "lipid A modification"),
              "DB"=c("COG","COG","COG","KEGG","KEGG","KEGG","KEGG","KEGG"),
              "Rel.Abun"=c(2.9,0.4,2.1,0.05,0.024,1.09,0.63, 0.043));

# plot
path.plot=ggplot(df,aes(x=Predictors, y=Rel.Abun,fill=DB))+
  geom_bar(stat="identity", 
           color="black")+ ##fcbba1
  facet_grid(~DB, scales="free",space = "free_x")+
  scale_y_continuous(limits = c(0, 5))+
  scale_fill_manual(values=c("palegreen","lightsalmon"))+
  labs(y="Relative Abundance (%)",
       x="")+
  theme_classic() + 
  theme(axis.title = element_text(size = 12, face="bold"), 
        axis.text = element_text(size = 12),
        axis.text.x=element_text(size=11, angle=90, vjust=0.5,face="bold"),
        legend.position="none");
plot(path.plot);

# save image 
ggsave(filename="07_vocKEGG.COG_pathways.png",
       device="png",path="./figures",
       plot=path.plot,
       width=4,
       height=10,
       units="in",
       dpi=500);



