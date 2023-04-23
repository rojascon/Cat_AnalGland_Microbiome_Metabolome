################################################################################
#
#                       Cat Anal Gland Microbiome 
#                      
#         Rojas et al 2023. Characterization of the microbiome and 
#         volatile compounds in anal gland secretions from domestic
#                 cats using metagenomics and metabolomics
#
#                       By: Connie Rojas
#                       Created: 16 August 2022
#                       Last updated: 20 April 2023
#
################################################################################

## CODE FOR: A) ploting MAG annotation rate
#            B) making Pie of MAG Genus-level abundances
#            C) making phylogenetic trees of each MAG alongside genomes from
#                             GTDB database
#            C) making barplot that shows bacterial species that were
#               recovered as cultured isolates, were present in the microbiome dataset
#               and had a reconstructed genome (MAGs)

source(file="scripts/00_configure.R"); #load necessary packages and specifications


################################################################################
#         1. Load MAG list, MAG abundances as calculated by CoverM software, 
#                               and sample metadata
################################################################################

# List of 89 quality MAGs
mag=read.csv("data/00_MAG_List.csv", header=T); 
mag=mag[mag$assembler=="metaSPADES",];

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
#                        2. Plot MAG Annotation Rate
#     # MAGs annotated at that classification level / total number of MAGs
################################################################################

# calculate annotation rates
categs=c("Phylum","Class","Order","Family","Genus","Species");
unclass=c(sum(is.na(mag$Phylum)),sum(is.na(mag$Class)), sum(is.na(mag$Order)),
          sum(is.na(mag$Family)),sum(is.na(mag$Genus)), sum(is.na(mag$Species)));

pers=100-((unclass/nrow(mag)) * 100);
andf=data.frame(categs,pers);
andf$categs=factor(andf$categs, levels=c("Phylum","Class","Order",
                                         "Family","Genus","Species"));
andf$categs=forcats::fct_rev(factor(andf$categs));

# plot!
magb=ggplot(data=andf, 
             aes(x=pers,y=categs))+
  geom_bar(stat = "identity",position="dodge",fill="#d4b9da",colour="black")+
  theme_bw()+ 
  labs(x = "Annotation Rate (%)",
       y = "")+
  theme(legend.position="none", 
        plot.title=element_text(size = 12, face="bold"),
        panel.background = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        axis.title.x = element_text(size=11, face="bold"),
        axis.text.x = element_text(size=12));plot(magb);

# save image
ggsave(filename="09_MAG_annotationRate.pdf",
       device="pdf",path="./figures",
       plot=magb,
       width=5,
       height=4,
       units="in",
       dpi=500);


################################################################################
#              3. Pie of MAG Family Abundances
################################################################################

# aggregate MAG abundances by Genus
rownames(mag)=mag$MAG.ID; mag$MAG.ID=NULL;
magab=merge(mag[,8:9],coverm,by=0);
magab=magab[,c(2,4:ncol(magab))];

magab$Genus[is.na(magab$Genus)]="non-Genus";
magab2=aggregate(.~Genus,magab, sum);

# calculate proportions
magab2[,-1] <- lapply(magab2[,-1], function(x) (x/sum(x))*100);
print(colSums(magab2[-1]));
magab2$AVG=rowMeans(magab2[,-1]);

# retain MAGs with mean abundances > 0.75
magab2=magab2[magab2$AVG>0.75,];
newrow=c(NA, 100-colSums(magab2[2:ncol(magab2)])); 
magab2=rbind(magab2, newrow); 
magab2[nrow(magab2),1]="Other";
magab2$sample="averages";

# select color-palette
piecol=c("maroon", "darkturquoise", "#225ea8", "gold1", "#FF7F00","yellow3",
         "brown", "gray70", "green4", "#FB9A99", "darkorange4", "deeppink1",
         "#6A3D9A", "steelblue4", "skyblue2", "yellow4", "palegreen2", 
         "black", "khaki2", "#CAB2D6", "#bf812d", "#FDBF6F");

# plot piechart
magp=ggplot(magab2, aes(x=sample, y=AVG, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=piecol)+
  theme_void()+
  labs(fill="",
       x="",
       title="MAG Genus Abundances")+
  theme(axis.text.x=element_blank(),
        text=element_text(size=10, face="bold"),
        legend.position = "bottom",
        legend.text = element_text(size=10),
        axis.title.y=element_text(size=12, face="bold"))+
  guides(fill = guide_legend(ncol=3));plot(magp);

# save image
ggsave(filename="09_MAG_Abun_Pie.pdf",
       device="pdf",path="./figures",
       plot=magp,
       width=5,
       height=4,
       units="in",
       dpi=500);


################################################################################
#         4.        Load GTDB phylogenetic tree and alignments
################################################################################

# load a formatted  list of the 89 quality MAGs so that it matches the formatting of the
# full GTDB taxonomic labels database (e.g. add c__ before each Bacterial Class label)
bin_names=read.csv("data/00_GTDB_MAG_tax_forTrees.csv", 
                   header=T);

# load full GTDB-r202 taxonomic labels downloaded from 
# https://data.gtdb.ecogenomic.org/releases/
load("data/00_GTDB_r202_taxonomy.Rdata");

# load alignments generated by GTDB when assigning taxonomy to MAGs
phylo.r202 <- read.newick("data/00_spades_gtdbtk.bac120.classify.tree");
phylo.r202$edge.length <- phylo.r202$edge.length*10;

# add your MAG taxonomic labels to the larger GTDB taxonomy
tree.meta.ms <- full_join(taxonomy.bac, bin_names);


################################################################################
#         5.  Subset GTDB phylogeny to only MAG of interest and replace 
#             default tree node labels with GTDB taxonomic labels
################################################################################

# select MAG of interest and how how tree levels to show
mag.tree <- tree_subset(phylo.r202, "metabat.19", levels_back = 4);

# subset larger GTDB taxonomy to only your labels of interest
mag.tree.taxaIDs <- mag.tree$tip.label;
tree.meta.ms2 <- tree.meta.ms %>% filter(label %in% mag.tree.taxaIDs);

# add taxonomic labels to your tree in preparation for plotting
my.tree <-full_join(mag.tree, tree.meta.ms2, by = "label");


################################################################################
#     6.  Plot phylogeny of MAG placed within the greater GTDB phylogeny
################################################################################

# plot tree
ptree=ggtree(my.tree, color = "black", size = 1, linetype = 1) + 
  geom_tiplab(aes(label = Taxonomy, color=Source),
              fontface = "bold.italic", size = 4, offset = 0.001) +
  xlim(0,6) + # (0,11) if detailed tree 
  scale_color_manual(values =c(study = "mediumorchid3", ref = "#000000")) +
  theme(legend.position = "none"); plot(ptree);

# save tree
ggsave(filename = 'figures/09_MAG19_phylogeny.png', 
       plot = last_plot(), 
       device = 'png', 
       width = 8, 
       height = 7, 
       dpi = 300);


################################################################################
#         7. Were there bacterial species that were recovered as cultured isolates,
#   were present in the microbiome dataset, AND had a reconstructed genome (MAG)?
################################################################################

# remove MAGs not classified to Species-level
msp=mag[complete.cases(mag$Species),];
msp=msp[,c(10,11,13)];
msp=msp[!duplicated(msp$final.name),];

# melt data frame for ggplot
mspp<-reshape2::melt(msp, id.vars="final.name",value.name = "appear");
mspp$variable=as.character(mspp$variable);
mspp$variable[mspp$variable=="sp_cultures"]="Cultured isolate";
mspp$variable[mspp$variable=="sp_braken"]="Microbiome data";

# plot
magplot=ggplot(data=mspp, aes(x=final.name, 
                            y=variable, fill=appear)) + 
  geom_tile() +
  facet_grid(~variable,scales="free_x",labeller = 
               labeller(variable = label_wrap_gen(width =12)))+ 
  scale_fill_manual(values=c("#f0f0f0","#7bccc4"), labels=c("N","Y")) +
  coord_flip() +
  labs(y="",
       x="",
       fill="")+
  theme_bw()+
  theme(legend.position="none",
        legend.title=element_text(size=13, face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10))+
  theme(strip.text.x = element_text(size = 11,face="bold")); 
plot(magplot); 

ggsave(filename="09_MAGs_Isolates_Microbiome.pdf",
       device="pdf",path="./figures",
       plot=magplot,
       width=5,
       height=7,
       units="in",
       dpi=500);