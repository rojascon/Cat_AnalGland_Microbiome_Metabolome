################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 15 Aug 2022
#                       Last updated: 21 April 2023
#
################################################################################

## CODE FOR: A) tidying data generated by Anvio after annotating contigs with
#           Cluster of Orthologous Genes (COG) & Kyoto Enc.Genes Genomes (KEGG)

#           B) appending gene counts to Anvio gene annotations
#           Gene counts were calculated using Salmon and are in TPM

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#             1. Load Anvio Gene annotations and Salmon Gene counts 
################################################################################
# load Salmon counts
load("data/00_Salmon_GeneTPM_SPADES.Rdata");
tpm$geneID=as.character(tpm$geneID);

# load Anvio Gene Annotations
load("data/00_Anvio_GeneAnnotations_SPADES.Rdata");

# only keep annotations that have statistically significant e.values
afun=afun[afun$e_value<=0.05,];


################################################################################
#             2. Tidy and subset KEGG Families (aka 'genes')          
################################################################################

# subset to only KEGG families
kegg=afun[afun$source=="KOfam",];row.names(kegg)=NULL;

# extract and tidy kegg annotations
kegg2=kegg[,c(1,4)];
kegg2b=as.data.frame(str_split_fixed(kegg2$genefunction, "EC",2));
kegg2$gf=kegg2b$V1;
kegg2$gf[kegg2$genefunction=="ECM component binding autotransporter adhesin"]=
  "ECM component binding autotransporter adhesin";

# append gene counts and aggregate by gene function
kc=merge(tpm, kegg2[,c(1,3)],by="geneID");
kc$geneID=NULL;
keggabun=aggregate(.~gf,kc, sum);  

#save for accessing later
save(keggabun, file="data/06_KOfam_abundances_spades.Rdata");


################################################################################
#             3. Tidy and subset KEGG Modules (aka 'pathways')          
################################################################################

# subset to only KEGG Modules
kegg=afun[afun$source=="KEGG_Module",];row.names(kegg)=NULL;

# extract and tidy kegg annotations
kegg2=kegg[,c(1,4)];
kegg2b=as.data.frame(str_split_fixed(kegg2$genefunction, "  ",2));
kegg2$gf=kegg2b$V1;

# append gene counts and aggregate by gene function
kc=merge(tpm, kegg2[,c(1,3)],by="geneID");
kc$geneID=NULL;
keggabun=aggregate(.~gf,kc, sum); 

#save for accessing later
save(keggabun, file="data/06_KOModule_abundances_spades.Rdata");


################################################################################
#             4. Tidy and subset COG functions (aka 'genes')          
################################################################################

# subset to only COG functions
cogc=afun[afun$source=="COG20_FUNCTION",];row.names(cogc)=NULL;

# extract and tidy COG annotations
cogc2=cogc[,c(1,4)];
cogc2b=as.data.frame(str_split_fixed(cogc2$genefunction, "   ",2));
cogc2$gf=cogc2b$V1;

# append gene counts
cc=merge(tpm, cogc2[,c(1,3)],by="geneID");
cc$geneID=NULL;
cogabun=aggregate(.~gf,cc, sum); # gf = gene function, cc= df of cog counts

# some manual editing 
cogabun$gf[cogabun$gf=="2  3"]="cNMP phosphodiesterase YmdB";
cogabun$gf[cogabun$gf=="3"]="ssDNA RNA exonuclease TatD";
cogabun$gf[cogabun$gf=="5  3"]="deoxyribonucleotidase  YorC";

#save for accessing later
save(cogabun, file="data/06_COGfun_spades.Rdata");


################################################################################
#             5. Tidy and subset COG Pathways     
################################################################################

# subset to only COG Pathways
cogc=afun[afun$source=="COG20_PATHWAY",];row.names(cogc)=NULL;

# extract and tidy COG annotations
cogc2=cogc[,c(1,4)];
cogc2b=as.data.frame(str_split_fixed(cogc2$genefunction, "   ",2));
cogc2$gf=cogc2b$V1;

# append gene counts
cc=merge(tpm, cogc2[,c(1,3)],by="geneID");
cc$geneID=NULL;
cogabun=aggregate(.~gf,cc, sum); 

#save for accessing later
save(cogabun, file="data/06_COGpath_spades.Rdata");
