#Btl Phylogenetic Tree Figures
#Morgan Carter 2023 R v 4.2.2
#Inputs are phylogenetic analyses as newick files
#Visualizes them with additional annotation based on modified info file from BtlRVDFinder


# Packages ----------------------------------------------------------------
setwd("") #Change for your use 
library(ggtree)
library(ggplot2)
library(dplyr)
library(treeio)
library(phytools)
library(ape)
library(ggpubr)
library(ggnewscale)
library(RColorBrewer)


# Inputs and Setup ------------------------------------------------------------------
tree_nwk_dis <- treeio::read.newick('MEGA/Intact_MEGA_Maximum_Likelihood.nwk', node.label='support')
df_meta<- read.csv("Btl_info_mod.csv", header=T)
#Info mod file is output from BtlRVDFinder with sequence source added by hand and any changes to protein names from manual curation

# Formatting Data ---------------------------------------------------------
df_meta$fastaName <- unlist(lapply(df_meta$fastaName, function(x) gsub(" ", "_", x)))
tree_nwk_dis@phylo[["tip.label"]] <- lapply(tree_nwk_dis@phylo[["tip.label"]], function(x) gsub(" ", "_", x))
tree_nwk_dis@phylo[["tip.label"]] <- unlist(lapply(tree_nwk_dis@phylo[["tip.label"]], function(x) gsub("'", "", x)))
setdiff(tree_nwk_dis@phylo[["tip.label"]], df_meta$fastaName) #check that they are the same


# Basic Tree Visualization ------------------------------------------------
tree_nwk_dis_visual <- ggtree(tree_nwk_dis, size=0.7)%<+%df_meta+# read in tree
  geom_tippoint(aes(color=Source),
                size=2)+
  scale_color_manual(
    name = "Source",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("NCBI", "ZygoLife", "PacBio"),                     # the different options in your variable
    values = c("#446600", "#9e0c42", "#5671d4")) +  
  new_scale_color()+
  theme(legend.title=element_text(face="bold"))+
  geom_tiplab(aes(label=Btl, color=Complete), size=3.5, align=TRUE, offset=0.01, show.legend = F)+
  scale_color_manual(values=c(Yes="black", No="grey"))+
  geom_nodelab(aes(label=round(support,2)*100, subset=support > 0.50), size=2.5, vjust=-0.5, hjust=1.4)
tree_nwk_dis_visual+theme_tree2()
tree1 <- tree_nwk_dis_visual + new_scale_fill()


# Adding Heatmaps ---------------------------------------------------------

#subsetting for heatmap data
df_meta_subset2 <- df_meta%>%select(Repeat_no)
rownames(df_meta_subset2) <- df_meta$fastaName
df_meta_subset0 <- df_meta%>%select(PseudoRepeat)
rownames(df_meta_subset0) <- df_meta$fastaName
df_meta_subset1 <- df_meta%>%select(FirstRepeat)
rownames(df_meta_subset1) <- df_meta$fastaName
df_meta_subset3 <- df_meta%>%select(LastRepeat)
rownames(df_meta_subset3) <- df_meta$fastaName
df_meta_subset4 <- df_meta%>%select(NLS)
rownames(df_meta_subset4) <- df_meta$fastaName
df_meta_subset5 <- df_meta%>%select(Species)
rownames(df_meta_subset5) <- df_meta$fastaName

#adding heatmaps
tree_heat0 <- 
  gheatmap(tree1, df_meta_subset0, offset=0.055, font.size = 4, width=0.04, color="grey", colnames=F)+
  scale_fill_manual(values=c("HY"="#ad224d", "QY"="#c76176", "HS"="#de96a1", "NY"="#f0cacf"), na.value = "white", name='Pseudo Repeat')+
  new_scale_fill()
tree_heat0
tree_heat1 <- 
  gheatmap(tree_heat0, df_meta_subset1, offset=0.07, font.size = 4, width=0.04, color="grey", colnames=F)+
  scale_fill_manual(values=c("YD"="#31410c", "HD"="#618218", "CD"="#b0c08b"), na.value = "white", name='First Repeat')+
  new_scale_fill()
tree_heat2 <- 
  gheatmap(tree_heat1, df_meta_subset2, offset=0.085, font.size = 4, width=0.04, hjust=0, color="grey", colnames=F)+
  scale_fill_gradient(low="#ffffff", high="#5671d4", name='# of Repeats')+
  ylim(NA, dim(df_meta)[1]+2)+
  new_scale_fill()
tree_heat3 <- 
  gheatmap(tree_heat2, df_meta_subset3, offset=0.1, font.size = 4, width=0.04, hjust=0, color="grey", colnames=F)+
  scale_fill_manual(values=c("N*"="#7b498b", "Other"="#a26eb3"), na.value = "white", name='Last Repeat')+
  new_scale_fill()
tree_heat4 <- 
  gheatmap(tree_heat3, df_meta_subset4, offset=0.115, font.size = 4, width=0.04, hjust=0, color="grey", colnames=F)+
  scale_fill_manual(values=c("RIRK"="#52906c", "QIRK"="#c4d9cc"), na.value = "white", name='NLS')+
  new_scale_fill()
tree_heat5 <-
  gheatmap(tree_heat4, df_meta_subset5, offset=0.13, font.size = 4, width=0.04, hjust=0, color="grey", colnames=F)+
  scale_fill_manual(values=c("M. rhizoxinica"="#f99708", "M. endofungorum"="#fbc06a", "M. sp."="#ffffff", "Unknown"="grey"), name='Species')+
  new_scale_fill()
tree_final <- tree_heat5+
  geom_cladelab(node=51, label="I", offset = 0.15, align = TRUE)+
  geom_cladelab(node=58, label="II", offset = 0.16, align = TRUE)+
  geom_cladelab(node=64, label="III", offset = 0.15, align = TRUE)+
  geom_cladelab(node=69, label="IV", offset = 0.16, align = TRUE)+
  geom_cladelab(node=78, label="V", offset = 0.15, align = TRUE)+
  geom_cladelab(node=83, label="VI", offset = 0.16, align = TRUE)
tree_final

#Export 10x10


# All Btl Tree ------------------------------------------------------------
# Repeats above but with the input of a larger tree

#load
tree_nwk_dis <- treeio::read.newick('MEGA/Btls_MEGA_Maximum_Likelihood_Condensed.nwk', node.label='support')

#add in data - make sure you have added a source and species column manually
df_meta<- read.csv("Btl_info_mod.csv", header=T)
df_meta$fastaName <- unlist(lapply(df_meta$fastaName, function(x) gsub(" ", "_", x)))
tree_nwk_dis@phylo[["tip.label"]] <- lapply(tree_nwk_dis@phylo[["tip.label"]], function(x) gsub(" ", "_", x))
tree_nwk_dis@phylo[["tip.label"]] <- unlist(lapply(tree_nwk_dis@phylo[["tip.label"]], function(x) gsub("'", "", x)))

setdiff(tree_nwk_dis@phylo[["tip.label"]], df_meta$fastaName) #check that they are the same

#visualize
tree_nwk_dis_visual <- ggtree(tree_nwk_dis, size=0.7)%<+%df_meta+# read in tree
  geom_tippoint(aes(color=Source),
                size=2)+
  scale_color_manual(
    name = "Source",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("NCBI", "ZygoLife", "PacBio"),                     # the different options in your variable
    values = c("#446600", "#9e0c42", "#5671d4")) +  
  new_scale_color()+
  #xlim(NA, 10)+
  theme(legend.title=element_text(face="bold"))+
  geom_tiplab(aes(label=Btl, color=Complete), size=3.5, align=TRUE, offset=0.01, show.legend = F)+
  scale_color_manual(values=c(Yes="black", No="grey"))+
  geom_nodelab(aes(label=round(support,2)*100, subset=support > 0.50), size=2.5, vjust=-0.5, hjust=1.4)
tree_nwk_dis_visual
tree1 <- tree_nwk_dis_visual + new_scale_fill()
#subsetting for heatmap data
df_meta_subset2 <- df_meta%>%select(Repeat_no)
rownames(df_meta_subset2) <- df_meta$fastaName
df_meta_subset0 <- df_meta%>%select(PseudoRepeat)
rownames(df_meta_subset0) <- df_meta$fastaName
df_meta_subset1 <- df_meta%>%select(FirstRepeat)
rownames(df_meta_subset1) <- df_meta$fastaName
df_meta_subset3 <- df_meta%>%select(LastRepeat)
rownames(df_meta_subset3) <- df_meta$fastaName
df_meta_subset4 <- df_meta%>%select(NLS)
rownames(df_meta_subset4) <- df_meta$fastaName
df_meta_subset5 <- df_meta%>%select(Species)
rownames(df_meta_subset5) <- df_meta$fastaName
#adding heatmaps
tree_heat0 <- 
  gheatmap(tree1, df_meta_subset0, offset=1, font.size = 4, width=0.02, color="grey", colnames=F)+
  scale_fill_manual(values=c("HY"="#ad224d", "QY"="#c76176", "HS"="#de96a1", "NY"="#f0cacf"), na.value = "white", name='Pseudo Repeat')+
  new_scale_fill()
tree_heat0
tree_heat1 <- 
  gheatmap(tree_heat0, df_meta_subset1, offset=1.15, font.size = 4, width=0.02, color="grey", colnames=F)+
  scale_fill_manual(values=c("YD"="#31410c", "HD"="#618218", "CD"="#b0c08b"), na.value = "white", name='First Repeat')+
  new_scale_fill()
tree_heat2 <- 
  gheatmap(tree_heat1, df_meta_subset2, offset=1.3, font.size = 4, width=0.02, hjust=0, color="grey", colnames=F)+
  scale_fill_gradient(low="#ffffff", high="#5671d4", name='# of Repeats')+
  ylim(NA, dim(df_meta)[1]+2)+
  new_scale_fill()
tree_heat3 <- 
  gheatmap(tree_heat2, df_meta_subset3, offset=1.45, font.size = 4, width=0.02, hjust=0, color="grey", colnames=F)+
  scale_fill_manual(values=c("N*"="#7b498b", "Other"="#a26eb3"), na.value = "white", name='Last Repeat')+
  new_scale_fill()
tree_heat4 <- 
  gheatmap(tree_heat3, df_meta_subset4, offset=1.60, font.size = 4, width=0.02, hjust=0, color="grey", colnames=F)+
  scale_fill_manual(values=c("RIRK"="#52906c", "QIRK"="#c4d9cc"), na.value = "white", name='NLS')+
  new_scale_fill()
tree_heat5 <-
  gheatmap(tree_heat4, df_meta_subset5, offset=1.75, font.size = 4, width=0.02, hjust=0, color="grey", colnames=F)+
  scale_fill_manual(values=c("M. rhizoxinica"="#f99708", "M. endofungorum"="#fbc06a", "M. sp."="#ffffff", "Unknown"="grey"), name='Species')+
  new_scale_fill()
tree_heat5


