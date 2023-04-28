################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 8 August 2022
#                       Last updated: 20 April 2023
#
################################################################################

## CODE FOR: A) determine whether microbiome beta-diversity varies with host
#               characteristics (age, diet, obesity, environment, periodontal disease)

#           B) construct PCoAs based on microbiota dissimilarity distances

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1. load taxa abundance table output by Kraken and sample metadata
################################################################################
load("data/01_bracken_Genus_nocontam.Rdata");
load("data/01_metadata_formatted.Rdata");

# subset to keep only anal sac samples 
meta2=meta[meta$swab_type=="Anal.Sac",];

# aggregate data by Genus
kspe=krak[,which(names(krak) 
                 %in% c(meta2$sampleID, "Genus"))];

colnames(kspe)[1]="taxa";
kspe=aggregate(.~taxa,kspe, sum);
rownames(kspe)=kspe$taxa; kspe$taxa=NULL;
kspe=kspe[rowSums(kspe)>0,];
kspe=kspe[,order(colnames(kspe))];


################################################################################
#             2. Transform data and then calculate 3 types of distances
################################################################################

# relative abundances for Bray-Curtis
bray<-apply(kspe, 2, function(i) (i/sum(i))*100);

# presence absence for Jaccard
jac=(kspe>0)*1;

# clr (center log ratio) for Aitchison distance ("compositions")
clra=data.frame(clr(kspe));

# generate distance matrices
ps1<- phyloseq(otu_table(bray, taxa_are_rows=TRUE));
bray.dist=phyloseq::distance(ps1, method="bray");

ps2<- phyloseq(otu_table(jac, taxa_are_rows=TRUE));
jac.dist=phyloseq::distance(ps2, method="jaccard",binary=T);

ps3<- phyloseq(otu_table(clra, taxa_are_rows=TRUE));
clr.dist=phyloseq::distance(ps3, method="euclidean");

# save the three distance matrices because will need them in future scripts
save("bray.dist", "jac.dist", "clr.dist", 
     file = "data/03_microbiome_taxa_distances.Rdata");


################################################################################
#             3. Generate Bray-Curtis distances and run PERMANOVAs 
#               age, diet, obesity, environment, periodontal disease)
################################################################################

# set variables for loop
mydist=list(jac.dist, bray.dist,clr.dist);
names=c("Jaccard","Bray-Curtis","Aitchison");
met=c("jaccard","bray","aitchison");

# run for loop to conduct the 3 PERMANOVA tests
for(i in 1:3)
{
  print(paste("PERMANOVA test, anal sac, N=23, using:", names[i]));
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
#             4. Plot PCoAs color-coded by age (yrs) or obesity (yes vs no)        
################################################################################

# calculate coordinates for PCoA
pcoa_dec=cmdscale(jac.dist, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");  
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sampleID");
pcoa_met=merge(pcoa,meta2,by="sampleID"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

# plot PCoA color-coded by age (yrs)
my_col=c("#253494","#7fcdbb","hotpink","gold");
pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=age_cat), 
             size = 2.7,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Age (yrs)")+
  scale_fill_manual(values=my_col)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa1);

# plot PCoA color-coded by obesity (Y vs N)
pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=obese), 
             size = 2.7,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Obesity")+
  scale_fill_manual(values=c("#c994c7","green"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa2);

# save images
age.bod=arrangeGrob(pcoa1,pcoa2, nrow=2);
ggsave(filename="03_pcoa_age.bodycondition.pdf",
       device="pdf",path="./figures",
       plot=age.bod,
       width=5,
       height=6.5,
       units="in",
       dpi=500);


