library(ape)
library(tidyverse)
library(ggtree)
library(ggpubr)

#Read in metadata
metadata <- read.csv("metadata.csv")


## GENOME SIZE FIGURES ##

#Add a column for cytometric completeness
metadata$cytometric.completeness <- metadata$assembly / metadata$flow.cytometry * 100

#Make dataframe of BUSCO and cytometric completeness
completeness.comp.df <- metadata %>%
  select(strain, taxon, assembly, flow.cytometry, busco.completeness, cytometric.completeness) %>%
  filter(!is.na(cytometric.completeness)) %>%
  gather(type, value, -strain, -taxon, -assembly, -flow.cytometry) %>%
  mutate(label=paste0('"("*italic("', taxon, '")*")"'))

#Format labels
completeness.comp.df$label <- sub(' sp\\.")\\*")', '")~"sp\\.)', completeness.comp.df$label)

#Barplot of completeness for each strain
gg.completeness <- ggplot(completeness.comp.df, aes(x=paste("IMI", strain), y=value, fill=type)) +
  geom_bar(stat="identity", 
           colour="black", 
           size=0.2,
           width=0.4,
           position="dodge") +
  geom_text(data=completeness.comp.df[completeness.comp.df$type == "busco.completeness",],
            aes(label=label, y=0),
            vjust=4,
            size=2,
            parse=T,
            show.legend=FALSE) +
  labs(y="Completeness (%)") +
  coord_cartesian(clip="off") +
  scale_y_continuous(limits=c(0, 100),
                     expand=c(0, 0)) +
  scale_fill_manual(values=c("white", "dimgrey"),
                    labels=c("Gene set (BUSCO)", "Cytometric genome size estimate"),
                    guide=guide_legend(title="Method of assessing assembly completeness",
                                       title.position="top",
                                       title.hjust=0.5)) +
  theme(legend.position="top",
        legend.key.size=unit(0.5,"line"),
        legend.title=element_text(size=9.5, face="bold"),
        legend.text=element_text(margin=margin(r=5)),
        legend.margin=margin(0,0,0,0),
        legend.box="vertical",
        axis.text.x=element_text(size=8, colour="black"),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_blank(),
        plot.margin=unit(c(1,1,5,1), "mm"),
        panel.grid.major.x=element_blank())

#Make dataframe of assembly and cytometry genome size
size.comp.df <- completeness.comp.df %>%
  select(strain, taxon, assembly, flow.cytometry, label) %>%
  distinct() %>%
  gather(type, value, -strain, -taxon, -label)

#Barplot of size for each strain
gg.size <- ggplot(size.comp.df, aes(x=paste("IMI", strain), y=value, fill=type)) +
  geom_bar(stat="identity", 
           colour="black", 
           size=0.2,
           width=0.4,
           position="dodge") +
  geom_text(data=size.comp.df[size.comp.df$type == "assembly",],
            aes(label=label, y=0),
            vjust=4,
            size=2,
            parse=T,
            show.legend=FALSE) +
  labs(y="Genome size (Mbp)") +
  coord_cartesian(clip="off") +
  scale_y_continuous(limits=c(0, 50),
                     expand=c(0, 0)) +
  scale_fill_manual(values=c("white", "dimgrey"),
                    labels=c("Total assembly size", "Cytometric genome size estimate"),
                    guide=guide_legend(title="Method of assessing genome size",
                                       title.position="top",
                                       title.hjust=0.5)) +
  theme(legend.position="top",
        legend.key.size=unit(0.5,"line"),
        legend.title=element_text(size=9.5, face="bold"),
        legend.text=element_text(margin=margin(r=5)),
        legend.margin=margin(0,0,0,0),
        legend.box="vertical",
        axis.text.x=element_text(size=8, colour="black"),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_blank(),
        plot.margin=unit(c(1,1,5,1), "mm"),
        panel.grid.major.x=element_blank())

#Write to file
#tiff(file=paste0("completeness_fig-", Sys.Date(), ".tiff"), height=5, width=6, units="in", res=300)
ggarrange(gg.size, gg.completeness, ncol=1, labels="AUTO", align="v")
#dev.off()


## GENETIC MARKERS FIGURE ##

#Make dataframe summarising markers for each lineage
markers.df <- metadata %>%
  select(lineage, acl1, Bml, cmdA, GAPDH, HIS3, mak.2, nik.1, PKC, LSU, RPB1, RPB2, TEF1, TUB2) %>%
  distinct() %>%
  gather(gene, copynum, -lineage) %>%
  mutate(gene=factor(gene)) %>%
  mutate(lineage=factor(lineage, levels=rev(levels(factor(lineage)))))

#Make vector to format lineage labels
label.face <- ifelse(is.na(str_extract(levels(markers.df$lineage), "aceae")), "italic", "plain")

#Grid plot of markers used for each lineage
gg.markers <- ggplot(markers.df, aes(x=gene, y=lineage, fill=copynum)) +
 geom_tile(colour="grey97", size=1) +
  scale_fill_manual(values=c("white", "dimgrey")) +
  scale_y_discrete(labels=str_to_sentence(levels(markers.df$lineage))) +
  scale_x_discrete(position="top",
                   labels=sub("\\.", "-", levels(markers.df$gene))) +
  theme_minimal() +
  theme(legend.position="none",
        legend.title=element_text(face="bold", hjust=1),
        axis.title=element_blank(),
        axis.text.x=element_text(size=6, colour="black"),
        axis.text.y=element_text(face=label.face, size=7, colour="black"))

#Write to file
#tiff(file=paste0("markers_fig-", Sys.Date(), ".tiff"), height=2, width=6, units="in", res=300)
gg.markers
#dev.off()


## T-BAS FIGURE ##

#Read in T-BAS tree
tbas <- read.tree("tbas_tree_220804.tre")

#Plot tree to get node numbers for clades containing own strains
#tiff(file=paste0("tbas-", Sys.Date(), ".tiff"), height=35, width=8, units="in", res=600)
ggtree(tbas, size=0.1) +
  geom_tiplab(size=0.8) +
  geom_text2(aes(subset=!isTip, label=node), size=0.8, hjust=-0.3)
#dev.off()

#Format labels
tbas$tip.label[which(tbas$tip.label == "355091")] <- "355091/355093"

#Make dataframe with nodes to plot lineages separately
tbas.df <- data.frame(node=c(3249, 3232, 3222, 3207, 2393, 2151, 2077, 2058, 2006),
           lineage=c("Didymellaceae", "Didymella", "Didymellaceae", "Didymosphaeriaceae", "Nectriaceae", "Colletotrichum", "Collariella", "Neurospora", "Gnomoniaceae"))

#Make vectors of strains
strains <- c("355080", "355082", "355084", "355091/355093", "356814", "356815", "359910", "360193", "360204", "364377", "366226", "366227", "366586", "367209")

#For each lineage...
for (i in 1:length(tbas.df$node)) {
  
  #Extract lineage clade
  tbas.tmp <- extract.clade(tbas, tbas.df$node[i])
  #Plot tree
  gg.tbas <- ggtree(tbas.tmp, size=0.3)
  
  #Create vector of fontface for new genomes in this study
  tip.face <- gg.tbas$data %>%
    filter(isTip == "TRUE") %>%
    pull(label)
  
  tip.face <- ifelse(tip.face %in% strains, "bold.italic", "italic")
  
  #Add formatted tip labels
  gg.tbas <- gg.tbas +
    xlim(0, 0.7) +
    geom_tiplab(size=2,
                font=tip.face,
                show.legend=FALSE) +
    ggtitle(label=tbas.df$lineage[i]) +
    scale_colour_manual(values=c("black", "red"))
  
  #Add lineage titles
  if (length(grep("aceae", tbas.df$lineage[i])) > 0) {
    
    gg.tbas <- gg.tbas +
      theme(plot.title=element_text(size=8, hjust=0.5, face="bold"))
    
  } else {
    
    gg.tbas <- gg.tbas +
      theme(plot.title=element_text(size=8, hjust=0.5, face="bold.italic"))
    
  }
  
  assign(paste0("gg.tbas.", i), gg.tbas)
  
}

#Write to file
#tiff(file=paste0("tbas_fig-", Sys.Date(), ".tiff"), height=5, width=6, units="in", res=600)
ggarrange(plotlist=mget(ls(pattern="gg.tbas.")))
#dev.off()


## TREE FIGURES ##

#Make list of outgroups for each tree
outgroups <- list(neocosmospora=c("Fusarium_staphyleae_strain_CBS_125482",
                                  "Fusarium_celtis-occidentalis_strain_CBS_125502",
                                  "Fusarium_cicatricum_strain_CBS_125550"),
                  colletotrichum=c("Monilochaetes_infuscans_culture-collection_CBS_869.96"),
                  collariella=c("Melanocarpus_albomyces_strain_CBS_638.94",
                                "Ovatospora_brasiliensis_strain_CBS_130174",
                                "Achaetomium_globosum_strain_CBS_332.67"),
                  neurospora=c("Sordaria_tomento-alba_CBS_260.78",
                               "Sordaria_fimicola_FGSC_2918",
                               "Sordaria_brevicollis_FGSC_1904",
                               "Pseudoneurospora_amorphoporcata_CBS_626.80",
                               "Sordaria_sclerogenia_FGSC_2741",
                               "Sordaria_macrospora_FGSC_4818"),
                  gnomoniopsis=c("Melanconis_alni_strain_AR_3500_voucher_BPI_748444",
                                 "Melanconis_marginalis_strain_AR_3442_voucher_BPI_748446"),
                  neocucurbitaria=c("Pyrenochaeta_acicola_CBS_812.95"),
                  neodidymelliopsis=c("Neoascochyta_desmazieri_strain_CBS_346.86",
                                      "Neoascochyta_sp._2_NV-2016",
                                      "Neoascochyta_europaea_strain_CBS_504.71"),
                  ascochyta=c("Phoma_herbarum_strain_CBS_615.75"),
                  didymosphaeriaceae=c("Massarinaceae_sp._DQD-2015a_voucher_MFLU_15-0057"),
                  didymella=c("Paraboeremia_putaminum_strain_CBS_130.69",
                              "Paraboeremia_adianticola_strain_CBS_187.83",
                              "Macroventuria_anomochaeta_strain_CBS_525.71",
                              "Macroventuria_wentii_strain_CBS_526.71",
                              "Paraboeremia_selaginellae_strain_CBS_122.93"))

#Make list of markers used for each tree
markers <- list(neocosmospora=c("acl1+cmdA+RPB1+\nRPB2+TEF1"),
                colletotrichum=c("GAPDH+HIS3+TUB2"),
                collariella=c("RPB2+TUB2"),
                neurospora=c("Bml+LSU+mak-2+\nnik-1+PKC+TEF1"),
                gnomoniopsis=c("TEF1+TUB2"),
                neocucurbitaria=c("RPB2+TEF1+TUB2"),
                neodidymelliopsis=c("RPB2+TUB2"),
                ascochyta=c("RPB2+TUB2"),
                didymosphaeriaceae=c("LSU+RPB2+\nTEF1+TUB2"),
                didymella=c("RPB2+TUB2"))

#Make list of axis limits for each tree
xlims <- list(neocosmospora=0.25,
              colletotrichum=0.7,
              collariella=0.4,
              neurospora=0.12,
              gnomoniopsis=0.5,
              neocucurbitaria=0.2,
              neodidymelliopsis=0.3,
              ascochyta=0.2,
              didymosphaeriaceae=0.6,
              didymella=0.25)

#Make list of tip label sizes for each tree
tip.sizes <- list(neocosmospora=3,
                  colletotrichum=2.3,
                  collariella=3,
                  neurospora=3,
                  gnomoniopsis=3,
                  neocucurbitaria=3,
                  neodidymelliopsis=3,
                  ascochyta=3,
                  didymosphaeriaceae=2.3,
                  didymella=3)

counter <- 0

#For each lineage...
for (i in c("didymella", "ascochyta", "neodidymelliopsis", "didymosphaeriaceae", "neocucurbitaria", "gnomoniopsis", "colletotrichum", "collariella", "neurospora", "neocosmospora")) {
  
  counter <- counter + 1
  
  #Read in tree
  tree <- read.tree(paste0("phylogenetics/raxmlng/", i, "/", i, "_concat.raxml.support"))
  #Root tree
  tree <- root(tree, outgroups[[i]], edgelabel=TRUE, resolve.root=TRUE)
  
  #If the lineage is one that contains an excessively long branch, truncate the branch
  if (i %in% c("neocosmospora", "neocucurbitaria", "neodidymelliopsis")) {
    
    shortened.edge <- tree$edge[which.max(tree$edge.length), 2]
    
    assign(paste0("shortened.edge.", i), shortened.edge)
    
    tree$edge.length[which.max(tree$edge.length)] <- tree$edge.length[which.max(tree$edge.length)] / 10
    
  }
  
  #Read in tree metadata
  metadata <- read.csv(paste0("phylogenetics/raxmlng/", i, "/metadata.csv"))
  
  #Plot tree
  gg.tree <- ggtree(tree, linetype=NA) %<+% metadata
  
  #Create vector of fontface for new assemblies in this study
  tip.face <- gg.tree$data %>%
    filter(isTip == "TRUE") %>%
    pull(own)
  
  tip.face <- ifelse(tip.face == "Y", "bold.italic", "italic")
  
  #Colour branches by support and scale x axis
  gg.tree <- gg.tree +
    geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
              show.legend=FALSE) +
    xlim(0, xlims[[i]]) +
    scale_colour_manual(values="darkgrey",
                        na.value="black")
  
  #If the lineage is one with species complexes to collapse...
  if (i %in% c("colletotrichum")) {

    #Make dataframe of species complex nodes
    clades.df <- data.frame(
      sc=unique(metadata %>%
                  filter(species.complex != "") %>%
                  pull(species.complex)),
      node=NA
    )

    #Find the most recent common ancestor for each clade
    for (j in 1:length(clades.df$sc)) {

      clades.df$node[j] <- MRCA(tree,
                                metadata$tip[metadata$species.complex == clades.df$sc[j]])

    }

    #Make dataframe of clades to collapse
    collapse.clades.df <- clades.df %>%
      filter(!sc %in% c("gloeosporioides", "acutatum")) %>%
      mutate(label=paste0("Colle. ", sc, " species complex"))

    for (j in 1:length(collapse.clades.df$sc)) {

      if (length(which(gg.tree$data$species.complex == collapse.clades.df$sc[j])) > 1) {
        
        #Collapse clades
        gg.tree <- ggtree::collapse(gg.tree,
                            node=collapse.clades.df$node[j])

      }

    }
    
    #Add clade labels
    gg.tree <- gg.tree +
      geom_cladelab(data=collapse.clades.df,
                    mapping=aes(node=node, label=label),
                    fontface="italic",
                    fontsize=tip.sizes[[i]],
                    offset.text=0.005) +
      geom_point2(aes(subset=node %in% collapse.clades.df$node),
                  shape=18, size=2)
    
    label.clades.df <- clades.df %>%
      filter(sc %in% c("gloeosporioides", "acutatum")) %>%
      mutate(label=paste0("Colle. ", sc, "\nspecies complex"))

    gg.tree <- gg.tree +
      geom_cladelab(data=label.clades.df,
                    mapping=aes(node=node, label=label),
                    fontface="italic",
                    fontsize=tip.sizes[[i]],
                    barsize=1,
                    offset=0.23,
                    offset.text=0.01)

  }
  
  #Add tip labels and a label indicating markers used
  gg.tree <- gg.tree +
    geom_tiplab(aes(label=name),
                size=tip.sizes[[i]],
                font=tip.face) +
    annotate(geom="text",
             label=LETTERS[counter],
              x=-Inf, y=Inf, hjust=0, vjust=1,
             size=5, fontface="bold") +
    annotate(geom="label",
             label=markers[[i]],
             x=max(na.omit(gg.tree$data$x))*0.17, y=max(na.omit(gg.tree$data$y))*0.92,
             size=3)
  
  assign(paste0("gg.tree.", i), gg.tree)

}

#Add labels indicating truncated branches
gg.tree.neocosmospora <- gg.tree.neocosmospora +
  geom_label2(aes(x=branch, subset=node == shortened.edge.neocosmospora),
              label="//",
              label.padding=unit(0, "pt"),
              label.size=0)

gg.tree.neocucurbitaria <- gg.tree.neocucurbitaria +
  geom_label2(aes(x=branch, subset=node == shortened.edge.neocucurbitaria),
              label="//",
              label.padding=unit(0, "pt"),
              label.size=0)

gg.tree.neodidymelliopsis <- gg.tree.neodidymelliopsis +
  geom_label2(aes(x=branch, subset=node == shortened.edge.neodidymelliopsis),
              label="//",
              label.padding=unit(0, "pt"),
              label.size=0)

#Write to file
#png(file="chapter4figure3.A.png", height=9, width=6, units="in", res=600)
#tiff(file=paste0("didymella_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.didymella
#dev.off()

#png(file="chapter4figure3.B.png", height=4, width=6, units="in", res=600)
#tiff(file=paste0("ascochyta_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.ascochyta
#dev.off()

#png(file="chapter4figure3.C.png", height=4.5, width=6, units="in", res=600)
#tiff(file=paste0("neodidymelliopsis_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.neodidymelliopsis
#dev.off()

#png(file="chapter4figure3.D.png", height=9, width=6, units="in", res=600)
#tiff(file=paste0("didymosphaeriaceae_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.didymosphaeriaceae
#dev.off()

#png(file="chapter4figure3.E.png", height=4.5, width=6, units="in", res=600)
#tiff(file=paste0("neocucurbitaria_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.neocucurbitaria
#dev.off()

#png(file="chapter4figure3.F.png", height=4.7, width=6, units="in", res=600)
#tiff(file=paste0("gnomoniopsis_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.gnomoniopsis
#dev.off()

#png(file="chapter4figure3.G.png", height=9, width=6, units="in", res=600)
#tiff(file=paste0("colletotrichum_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.colletotrichum
#dev.off()

#png(file="chapter4figure3.H.png", height=4.5, width=6, units="in", res=600)
#tiff(file=paste0("collariella_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.collariella
#dev.off()

#png(file="chapter4figure3.I.png", height=6.5, width=6, units="in", res=600)
#tiff(file=paste0("neurospora_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.neurospora
#dev.off()

#png(file="chapter4figure3.J.png", height=9, width=6, units="in", res=600)
#tiff(file=paste0("neocosmospora_fig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600)
gg.tree.neocosmospora
#dev.off()

#Read in supplementary LSU tree
supp.tree <- read.tree("phylogenetics/raxmlng/didymosphaeriaceae/didymosphaeriaceae_LSU.raxml.support")
#Root
supp.tree <- root(supp.tree, outgroups$didymosphaeriaceae, edgelabel=TRUE, resolve.root=TRUE)
#Truncate branch
shortened.edge <- supp.tree$edge[which.max(supp.tree$edge.length), 2]
supp.tree$edge.length[which.max(supp.tree$edge.length)] <- supp.tree$edge.length[which.max(supp.tree$edge.length)] / 3
#Read in metadata
metadata <- read.csv("phylogenetics/raxmlng/didymosphaeriaceae/metadata.csv")

#Plot tree
gg.supp.tree <- ggtree(supp.tree, linetype=NA) %<+% metadata

#Create vector of fontface for new genomes in this study
tip.face <- gg.supp.tree$data %>%
  filter(isTip == "TRUE") %>%
  pull(own)

tip.face <- ifelse(tip.face == "Y", "bold.italic", "italic")

#Plot tree
gg.supp.tree <- gg.supp.tree +
  geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  xlim(0, xlims[[i]]) +
  scale_colour_manual(values="darkgrey",
                      na.value="black") +
  geom_label2(aes(x=branch, subset=node == shortened.edge),
              label="//",
              size=3,
              label.padding=unit(0, "pt"),
              label.size=0) +
  xlim(0, 0.15) +
  geom_tiplab(aes(label=name),
              size=2.2,
              font=tip.face)

#Write to file
#tiff(file=paste0("didymosphaeriaceae_suppfig-", Sys.Date(), ".tiff"), height=9, width=6, units="in", res=600, compression="lzw")
gg.supp.tree
#dev.off()