################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 26 Sept 2022
#                       Last updated: 20 April 2023
#
################################################################################

## CODE FOR: determining whether Kraken2 microbiome profiles are correlated with
#                   metabolome profiles extracted during solid-phase or
#                     after liquid derivatization using GC-MS

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1.  Load metabolite abundance table, and sample metadata
################################################################################

voc=read.csv("data/00_volatiles_solid_avg.csv", header=T); # solid-phase
lvoc=read.csv("data/00_volatiles_liquid_avg.csv", header=T); # liquid deriv.
load("data/01_metadata_formatted.Rdata");

# retain anal gland samples
meta2=meta[meta$swab_type=="Anal.Sac",];


################################################################################
#             2.  normalize, log-transform, and scale volatile data
#                           OR
#                 simply do presence/absence
################################################################################

# pick one (either solid-phase or liquid-phase);
mydf=voc[,c("compound",meta2$sampleID)];  
mydf=lvoc[,c("compound",meta2$sampleID)];  

# NORMALIZE with total sum scaling
mydf[,-1] <- lapply(mydf[,-1], function(x) (x/sum(x))*100);
print(colSums(mydf[-1]));

# LOG-10 transform
rownames(mydf)=mydf$compound; mydf$compound=NULL;
mydf=mydf+0.00000001;
mydf=log10(mydf);

# pareto scaling with mean center
mydf=data.frame(t(mydf));
mydf2<- data.frame(apply(mydf, 2, function(x) x - mean(x)));
mydf3<- data.frame(apply(mydf2, 2, function(x) x/sqrt(sd(x))));

# NOW PRESENCE/ABSENCE
# pick one (either solid-phase or liquid-phase);
jdf=voc[,c("compound",meta2$sampleID)];
jdf=lvoc[,c("compound",meta2$sampleID)];

# presence absence for jaccard 
rownames(jdf)=jdf$compound; jdf$compound=NULL;
jdf2=(jdf>0)*1;
jdf2=jdf2[,order(colnames(jdf2))];


################################################################################
#             3.  calculate distance matrices for volatile data
#                            and microbiome data
################################################################################

# volatile data euclidean distances
mydf3=mydf3[order(rownames(mydf3)),];
veuc.dist=stats::dist(mydf3, method = "euclidean");

# volatile data jaccard distances
ps<- phyloseq(otu_table(jdf2, taxa_are_rows=TRUE));
vjac.dist=phyloseq::distance(ps, method="jaccard",binary=T);

# save the distances because will need them in other scripts
save("vjac.dist", "veuc.dist", 
     file = "data/04_metabolite_solid_distances.Rdata");

# repeat code in section 2 and 3 for volatiles obtained after liquid-derivatization
save("vjac.dist", "veuc.dist", 
     file = "data/04_metabolite_liquid_distances.Rdata");

# load microbiome distances you already calculated in script 03_kraken_betadiversity.R
load("data/03_microbiome_taxa_distances.Rdata");


################################################################################
#             4.  Run the correlations
################################################################################

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


################################################################################
#            4.  plot metabolome principal coordinates ~ microbiome principal
#                     coordinates to illustrate relationship
################################################################################
load("data/04_metabolite_solid_distances.Rdata");

# generate volatile principal coordinates
pcoa_dec2=cmdscale(veuc.dist, eig=TRUE);  
pcoa2=as.data.frame(pcoa_dec2$points);
colnames(pcoa2)=c("Axis1","Axis2");

# generate microbiome principal coordinates
pcoa_dec=cmdscale(clr.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");

# plot volatile PC1 ~ microbiome PC1
cdf=data.frame(pcoa$Axis1, y=pcoa2$Axis1);
colnames(cdf)=c("Axis1","PC1");

vmp1=ggplot(cdf, aes(x=Axis1, y=PC1))+
  geom_point(color="black", size=2)+
  geom_smooth(method=lm, color="grey80", se=FALSE)+
  theme_classic()+labs(title = "")+
  labs(y="Metabolite PC1",
       x="Microbiome PC1")+
  theme(legend.title=element_blank(),text = element_text(size=12),
        legend.position="right",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold")); plot(vmp1);

# plot volatile PC2 ~ microbiome PC2
cdf=data.frame(pcoa$Axis2, y=pcoa2$Axis2);
colnames(cdf)=c("Axis2","PC2");

vmp2=ggplot(cdf, aes(x=Axis2, y=PC2))+
  geom_point(color="black", size=2)+
  geom_smooth(method=lm, color="grey80", se=FALSE)+
  theme_classic()+labs(title = "")+
  labs(y="Metabolite PC2",
       x="Microbiome PC2")+
  theme(legend.title=element_blank(),text = element_text(size=12),
        legend.position="right",
        plot.title = element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12, face="bold"),
        axis.title.y=element_text(size=12, face="bold")); plot(vmp2);

# save images
vmplot=arrangeGrob(vmp1,vmp2, nrow=1);
ggsave(filename="04_voc_microbiome_scatter.pdf",
       device="pdf",path="./figures",
       plot=vmplot,
       width=8,
       height=3.5,
       units="in",
       dpi=500);