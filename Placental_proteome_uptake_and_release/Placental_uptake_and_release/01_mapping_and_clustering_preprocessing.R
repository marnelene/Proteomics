###
#
# author: Maren-Helene Høie Degnes
# date: 04.11.21
# purpose: cluster vessels across all their proteins to identify 
#          subgroupings and map proteins to reactome categories to 
#          get an overview of the sub set of proteins we have 
#          with our 5000 proteins
#
###



library(tidyverse)
library(pls)
library(dendextend)
library(RColorBrewer)




# Load data
f.path <- "/data/lgRFU.RData"
load(file=f.path)






############### ================ cluster all log2 RFU samples from 4V cohort ============= #############
# sub set data
samprot.4v <- samprot[which(!is.na(samp$SampleGroup) & samp$Control==1),]
samp.4v <- samp[which(!is.na(samp$SampleGroup) & samp$Control==1),]

# cluster sub data 
crs <- cor(samprot.4v %>% t(), method = "pearson")
distance <- as.dist(1-(crs))
hclustering <- distance %>% hclust(.)
hclustering$labels <- ""                      # too many names to plot so we skip names of leaves
dend <- hclustering %>% as.dendrogram()
## from documentation for colored_bars
# As long as the sort_by_labels_order parameter is TRUE (default), the colors
# vector/matrix should be provided in the order of the original data order
# (and it will be re-ordered automatically to the order of the dendrogram)
colorvector <- c("red","green3","black","white","orange")
vesselnames <- c("Radial artery",
                 "Antecubital vein","Uterine vein",
                 "Umbilical artery","Umbilical vein")
grcolrs <- colorvector[samp.4v$SampleGroup %>% as.factor() %>% as.numeric()]


crs.mat <- crs[which(samp.4v$SampleGroup %in% c("PVE","UE","PAE")),which(samp.4v$SampleGroup %in% c("PVE","UE","PAE"))]
crs.mat <- crs.mat[upper.tri(crs.mat)]
mean(crs.mat); sd(crs.mat)

crs.fet <- crs[which(samp.4v$SampleGroup %in% c("AUE","VUE")),which(samp.4v$SampleGroup %in% c("AUE","VUE"))]
crs.fet <- crs.fet[upper.tri(crs.fet)]
mean(crs.fet);  sd(crs.fet)

crs.mat.fet <- crs[which(samp.4v$SampleGroup %in% c("PVE","UE","PAE")),which(samp.4v$SampleGroup %in% c("AUE","VUE"))]
crs.mat.fet <- crs.mat.fet[upper.tri(crs.mat.fet)]
mean(crs.mat.fet);  sd(crs.mat.fet)


tiff("SupplementaryFigure2_dend4V.tiff",
     units="px", width=2250, height=900, res=300, compression = 'lzw')
par(mfrow=c(1,1),mar=c(3,4,0,3))
plot(dend,yaxt="n",ylab="Pearson correlations")
axis(2,at=seq(0,0.3,by=0.05),labels=seq(1,0.7,by=-0.05))
colored_bars(colors = cbind(samp.4v$ID %>% as.factor() %>% as.numeric(),grcolrs),
             dend = dend,
             rowLabels = c("Patients","Vessel"),
             sort_by_labels_order = TRUE)

legend("topright",                    # With a specific distance from the plot margine
       inset = c(-0.00001, 0),
       title = "Vessel",        
       legend = vesselnames,
       fill = colorvector[c(2:4,1,5)],
       bty = "n",                             # No border around the legend
       title.adj = c(0.1))                    # Adjust the position of the legend title


dev.off()










################## ============== mapping of proteins to reactome categories =============== ##############

# Based on information found here https://reactome.org/download-data?id=61&ml=1 
# I can detect which pathways are parent and child in the uniprot-mapping data frame below!
reac.rela <- read_tsv("https://reactome.org/download/current/ReactomePathwaysRelation.txt", col_names = F)
colnames(reac.rela) <- c("parent", "child")

# Read in mappings between uniprot codes and reactome pathways 
reac.all.levels <- read_tsv("https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt", col_names = F)
colnames(reac.all.levels) <- c("Uniprot", "ReacPathwayID", "PathwayURL", "PathwayName", "Source", "Organism")
reac.alllv.hs <- reac.all.levels[which(reac.all.levels$Organism=="Homo sapiens"),]


# Create a list with protein names of proteins belonging to each pathway in each listbag
path <- c()
unique.pathways <- unique(reac.alllv.hs$ReacPathwayID)
genesets <- vector(mode = "list", length = length(unique.pathways))
pathnames <- c()
names(genesets) <- unique.pathways
for (pathnr in 1:length(unique.pathways)){
  protein.idx <- which(reac.alllv.hs$ReacPathwayID==unique.pathways[pathnr])
  genesets[[pathnr]] <- reac.alllv.hs$Uniprot[protein.idx]
  pathnames <- c(pathnames, unique(reac.alllv.hs$PathwayName[protein.idx]))
}



# Load pathway id dataframe of the 27 reactome pathway categories and prepare data frames
## Reactome categories and their ids
top.ids.df <- data.frame('Metabolism' = 1430728,
                         'Autophagy' = 9612973,
                         'Cell Cycle' = 1640170,
                         'Cell-Cell communication' = 1500931,
                         'Cellular responses to external stimuli' = 8953897,
                         'Chromatin organization' = 4839726,
                         'Circadian Clock' = 400253,
                         'Developmental Biology' = 1266738,
                         'Digestion and absorption' = 8963743,
                         'Disease' = 1643685,
                         'DNA Repair' = 73894,
                         'DNA Replication' = 69306,
                         'Extracellular matrix organization' = 1474244,
                         'Gene expression (Transcription)' = 74160,
                         'Hemostasis' = 109582,
                         'Immune System' = 168256,
                         'Metabolism of proteins' = 392499,
                         'Metabolism of RNA' = 8953854,
                         'Muscle contraction' = 397014,
                         'Neuronal System' = 112316,
                         'Organelle biogenesis and maintenance' = 1852241,
                         'Programmed Cell Death' = 5357801,
                         'Protein localization' = 9609507,
                         'Reproduction' = 1474165,
                         'Signal Transduction' = 162582,
                         'Transport of small molecules' = 382551,
                         'Vesicle-mediated transport' = 5653656)
top.ids <- rbind(top.ids.df, paste0("R-HSA-", top.ids.df[1,]), rep(NA, 27), rep(NA, 27),
                 rep(NA, 27), rep(NA, 27), rep(NA, 27),rep(NA, 27),rep(NA, 27),rep(NA, 27))
rownames(top.ids) <- c("ReactomeID", "ReactomeNameID", "PathwaySize","CountOfSomasInCategory", "CountOfMsec", 
                       "CountOfMupt", "CountOfFsec","CountOfFupt", "CountsOfArmsec","CountsOfArmupt")


## Create mapping function
# modified 07.12.20
# input: vector of uniprot ids
# output: vector of length n.paths with how many of input uniprots belonging to each pathway
n.paths <- 27
MapProteinSets <- function(uniprotvector){
  cat(length(uniprotvector), "\n")
  pathpresence <- c()
  pathsizes <- c()
  
  for (pathwaycolumn in 1:n.paths){
    pathwayid <- top.ids[2,pathwaycolumn]
    genesetidx <- which(names(genesets)==pathwayid)
    if(length(genesetidx)>1){cat("More than one pathway id in gene sets matched with current pathway id")
    } else if (length(genesetidx)<1){cat("No pathway id in gene sets mapped to top ids!")
    } else{
      pathpresence <- c(pathpresence, (which(genesets[[genesetidx]] %in% uniprotvector) %>% length()))
      pathsizes <- c(pathsizes, length(genesets[[genesetidx]]))
    }
  }
  presencedf <- data.frame(pathpresence=pathpresence,
                           pathsizes=pathsizes)
  return(presencedf)
}

# load and define uniprot vector containing all uniprot codes in the study
load("data/DupMultFixed.RData")  #load UniprotFixVector in same order as samprot. With NA for duplicated proteins. 
uniprotvec <- UniprotFixVector[-which(is.na(UniprotFixVector))]  # uniprotfixvector indicates duplicates uniprot codes with NA
presencemx <- MapProteinSets(uniprotvector=uniprotvec) %>% as.matrix() %>% t()   # formatting df to go into the barplot function properly
colnames(presencemx) <- top.ids.df %>% names()
tiff('SupplementaryFigure3_SOMAmap.tiff',
     units="px", width=2250, height=1700, res=300, compression = 'lzw',family = "sans")
par(mar=c(4.1,17,0,1))
barplot(presencemx, 
        horiz = T, 
        las=2,
        xlab="Number of proteins in Reactome Category")
legend("topright",
       fill=c("grey","gray35"),
       legend=c("Proteins associated",
                "Our proteins mapped"),
       bty="n");  dev.off()




