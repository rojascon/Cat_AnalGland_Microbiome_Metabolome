################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 29 Nov 2022
#                       Last updated: 21 April 2023
#
################################################################################

## CODE FOR: determining whether the relative abundaces of bacterial species
#                 vary with the relative abundances of identified and
#                     unidentified metabolites measured using GC-MS solid-phase microextraction

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1.  Load metabolite abundance table, microbiome species abundances, 
#                                 and sample metadata
################################################################################

voc=read.csv("data/00_volatiles_solid_avg.csv", header=T);
load("data/01_bracken_Genus_nocontam.Rdata");
load("data/01_metadata_formatted.Rdata");

# retain anal gland samples
meta2=meta[meta$swab_type=="Anal.Sac",];


################################################################################
#             2.  Reshape metabolite data in preparation for correlations
################################################################################

# calculate relative abundances and log10-transform
mydf=voc[,c("compound",meta2$sampleID)];  
mydf[,-1] <- lapply(mydf[,-1], function(x) (x/sum(x))*100);
rownames(mydf)=mydf$compound; mydf$compound=NULL;
mydf=mydf+0.00000001;
mydf=log10(mydf);

# pareto scaling
mydf=data.frame(t(mydf));
mydf2<- data.frame(apply(mydf, 2, function(x) x - mean(x)));
mydf3<- data.frame(apply(mydf2, 2, function(x) x/sqrt(sd(x))));
df4=data.frame(t(mydf3));

# keep all identified metabolites and only keep unidentified metabolites with abundances > 0
Ad=df4[1:37,];  # 
Bd=df4[38:nrow(df4),];
Bd=Bd[rowMeans(Bd)>0,]; # 
df5=rbind(Ad,Bd);
vocdata=df5[,order(colnames(df5))];


################################################################################
#             3.  Reshape microbiome data in preparation for correlations
################################################################################
# calculate relative abundances
kspe=krak[,which(names(krak) 
                 %in% c(meta2$sampleID, "Species"))]; 
colnames(kspe)[1]="taxa";
kspe=aggregate(.~taxa,kspe, sum);
rownames(kspe)=kspe$taxa; kspe$taxa=NULL;
kspe=kspe[rowSums(kspe)>0,];

# keep taxa with prevalence > 90%
gpad=(kspe>0)*1; 
gpad=as.data.frame(gpad);
corenum=round(0.90*ncol(gpad));  
gcor=gpad[rowSums(gpad)>corenum,];

# and keep taxa with abundance > 0.2
gd2=apply(kspe, 2, function(i) (i/sum(i))*100);
gd2=gd2[rownames(gd2) %in% rownames(gcor),];
gd2=gd2[rowMeans(gd2) > 0.2,]; 
gd2=data.frame(gd2);
bacdata=gd2[,order(names(gd2))];


################################################################################
#             4.  Correlate the 40 bacterial species with 37 identified metabolites
#                         and 193 unidentified metabolites
################################################################################
# get the names of the 40 bacterial species
rownames(bacdata);

# manually run this for each bacterial species (1-40)
# have to change the number inside as.numeric(bacdata[])
# in this example we have 40 for the 40th bacterial species
gres2=c();
for(i in 1:nrow(vocdata))
{
  myt=cor.test(as.numeric(vocdata[i,]),as.numeric(bacdata[40,]), method="spearman");   
  estim= myt$estimate
  pval= myt$p.value;
  tempdf=data.frame(rownames(vocdata)[i],estim, pval);
  gres2=rbind(gres2,tempdf);
};

gres2$padj=p.adjust(gres2$pval, method="fdr", n=length(gres2$pval));

# save each output manually
write.csv(gres2,"Trueperella pyogenes.csv",row.names=F);


################################################################################
#             5.  Only keep significant correlations (p>0.05)
################################################################################

# I compiled the results from the 40 bacterial species into one Excel file
# one file has the rho values, the other has the adjusted p.values
rhos=read.csv("data/05_rho_bacteria_metabolites.csv",row.names = 1, header=T);
pvals=read.csv("data/05_padj_bacteria_metabolites.csv",row.names = 1, header=T);

# remove non-sig correlations
rhos[pvals>0.05]="NA";
pvals[pvals>0.05]="NA";
View(rhos); View(pvals);

# save the significant correlations to their own files
# write.csv(rhos,file="rhos_significant.csv",row.names = T);
# write.csv(pvals,file="pvals_significant.csv",row.names = T);
