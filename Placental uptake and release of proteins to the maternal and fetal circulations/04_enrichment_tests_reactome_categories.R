###
#
# author: Maren-Helene Høie Degnes
# date: 04.11.21
# purpose: perform enrichment tests for the six classes of proteins:
#          - proteins released from placenta to maternal circulation
#          - proteins taken up from maternal circulation to placenta
#          - proteins released from placenta to fetal circulation
#          - proteins taken up from fetal circulation to placenta
#          - proteins upregulated released from placenta to maternal circulation
#          - proteins specifically released from placenta to maternal circulation
#
###




library(readr)
library(tidyverse)
library(ComplexHeatmap)


# Perform tests for specifically and/or upregulated released proteins?
placenta_specific <- T




################ ======================== Generate Genesets List ~ gmt obj   ====================== ################### 
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

plot(sapply(genesets, length), xlab="pathway number", ylab="Number of uniprot codes in pathway")





############# =================== Save Reactome categories names and ids ============ ###############
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





#################### ================= Map function ========== ################
MapProteinSets <- function(uniprotvector){
  cat(length(uniprotvector), "\n")
  pathpresence <- c()
  
  for (pathwaycolumn in 1:n.paths){
    pathwayid <- top.ids[2,pathwaycolumn]
    genesetidx <- which(names(genesets)==pathwayid)
    if(length(genesetidx)>1){cat("More than one pathway id in gene sets matched with current pathway id")
    } else if (length(genesetidx)<1){cat("No pathway id in gene sets mapped to top ids!")
    } else{pathpresence <- c(pathpresence, (which(genesets[[genesetidx]] %in% uniprotvector) %>% length()))
    }
  }
  return(pathpresence)
}






################### ================ Load data and test results for classes of proteins =============== ################
load("data/lgRFU.RData")
load("data/MTT_FBmedians.RData")
uniprotvec = prot$UniProt                   # MatNormTest$Uniprot contain all unique, single uniprot codes
allSOMAprotPathCount <- MapProteinSets(uniprotvector=uniprotvec)

# create
if(placenta_specific){
  # This is really specifically released proteins to maternal side, 
  # but for simplicity I give the same name as for if not placenta_specific
  matsec <- which(MatNormTest$`adjusted p`<0.05
                    & MatNormTest$`t statistic`>0
                    & ((ArmNormTest$`p value`>0.05)
                       | ArmNormTest$`t statistic`<0))
  # This is really upregulated released proteins to maternal side, 
  # but for simplicity I give the same name as for if not placenta_specific
  matupt <- which(MatNormTest$`adjusted p`<0.05
                   & MatNormTest$`t statistic`>0
                   & MatNormTest$PLARMpadj<0.05
                   & MatNormTest$PLARMt>0)
} else {
  matsec <- which(MatNormTest$`adjusted p`<0.05 & MatNormTest$`t statistic`>=0)
  matupt <- which(MatNormTest$`adjusted p`<0.05 & MatNormTest$`t statistic`<0)
  }



# Categorize results and save indeses for significant and categorized proteins
fetsec <- which(FetNormTest$`adjusted p`<0.05 & FetNormTest$`t statistic`>=0)
fetupt <- which(FetNormTest$`adjusted p`<0.05 & FetNormTest$`t statistic`<0)
armsec <- which(ArmNormTest$`adjusted p`<0.05 & ArmNormTest$`t statistic`>=0)
armupt <- which(ArmNormTest$`adjusted p`<0.05 & ArmNormTest$`t statistic`<0)






################ ================== Perform mapping prior to enrichment tests ======== ###########
n.paths <- ncol(top.ids.df)

# Load pathway id dataframe of the 27 reactome pathway categories and prepare data frames
load("K:/Sensitivt/Forskning01/2012-5678_Placenta/Maren-Helene/Rscripts/SOMA20/data/Reactop.RData")
top.ids <- rbind(top.ids.df, 
                 paste0("R-HSA-", top.ids.df[1,]), 
                 rep(NA, n.paths), rep(NA, n.paths), rep(NA, n.paths), rep(NA, n.paths), 
                 rep(NA, n.paths),rep(NA, n.paths),rep(NA, n.paths),rep(NA, n.paths))
rownames(top.ids) <- c("ReactomeID", "ReactomeNameID", "PathwaySize","CountOfSomasInCategory", "CountOfMsec", 
                       "CountOfMupt", "CountOfFsec","CountOfFupt", "CountsOfArmsec","CountsOfArmupt")
top.ids.coverage <- rbind(top.ids.df, 
                          paste0("R-HSA-", top.ids.df[1,]), 
                          rep(NA, n.paths),  rep(NA, n.paths), rep(NA, n.paths), 
                          rep(NA, n.paths),rep(NA, n.paths),rep(NA, n.paths), rep(NA, n.paths))
rownames(top.ids.coverage) <- c("ReactomeID", "ReactomeNameID","CovOfSomasInCategory","CovOfMsec", 
                                "CovOfMupt", "CovOfFsec","CovOfFupt","CovOfArmsec","CovOfArmupt")

# Save count of significant proteins in each category
top.ids$TestCount <- c(NA, NA, NA, NA, length(matsec), length(matupt), length(fetsec), length(fetupt), length(armsec), length(armupt))

# Save count of significant and categorized proteins present in protein sets for each pathway
top.ids[4,-28] <- MapProteinSets(uniprotvector = prot$UniProt)
top.ids[5,-28] <- MapProteinSets(uniprotvector = prot$UniProt[matsec])
top.ids[6,-28] <- MapProteinSets(uniprotvector = prot$UniProt[matupt])
top.ids[7,-28] <- MapProteinSets(uniprotvector = prot$UniProt[fetsec])
top.ids[8,-28] <- MapProteinSets(uniprotvector = prot$UniProt[fetupt])
top.ids[9,-28] <- MapProteinSets(uniprotvector = prot$UniProt[armsec])
top.ids[10,-28] <- MapProteinSets(uniprotvector = prot$UniProt[armupt])

# Save count of significant and categorized proteins present in protein sets for each pathway divided by amount of soma proteins in that pathway category
top.ids.coverage[3,] <- (MapProteinSets(uniprotvector = prot$UniProt)  / allSOMAprotPathCount) * 100
top.ids.coverage[4,] <- (MapProteinSets(uniprotvector = prot$UniProt[matsec])  / allSOMAprotPathCount) * 100
top.ids.coverage[5,] <- (MapProteinSets(uniprotvector = prot$UniProt[matupt])  / allSOMAprotPathCount) * 100
top.ids.coverage[6,] <- (MapProteinSets(uniprotvector = prot$UniProt[fetsec])  / allSOMAprotPathCount) * 100
top.ids.coverage[7,] <- (MapProteinSets(uniprotvector = prot$UniProt[fetupt])  / allSOMAprotPathCount) * 100
top.ids.coverage[8,] <- (MapProteinSets(uniprotvector = prot$UniProt[armsec])  / allSOMAprotPathCount) * 100
top.ids.coverage[9,] <- (MapProteinSets(uniprotvector = prot$UniProt[armupt])  / allSOMAprotPathCount) * 100

# save size of pathway (length of protein set independent on our soma proteins selsection) for all pathway categories
for (column in 1:n.paths){top.ids[3,column] <-  genesets[[which(names(genesets)==top.ids[2,column])]] %>% length()}   

# Create matrix to fill with p values
pvaluestab <-  matrix(NA,nrow=n.paths,ncol = 6)
rownames(pvaluestab) <- colnames(top.ids.df)
colnames(pvaluestab) <- paste0("p",c("matsec","matupt","fetsec","fetupt","armsec","armupt"))

# Create matrix to fill with odds ratios 
direction_of_fisher_tab <-  matrix(NA,nrow=n.paths,ncol = 6)
rownames(direction_of_fisher_tab) <- colnames(top.ids.df)
colnames(direction_of_fisher_tab) <- paste0("p",c("matsec","matupt","fetsec","fetupt","armsec","armupt"))





################## ================ Perform fishers's exact test ================= ############
# for each multiple testing results category for each pathway category
TotalNumberOfSOMAProteins <- length(uniprotvec)
NumberOfMultestResCategories <- ncol(pvaluestab)
for (MultestResCati in 1:NumberOfMultestResCategories){
  NumberOfProteinsSignInResCat <- top.ids$TestCount[MultestResCati+4]               # number of proteins in res class
  for (pathi in 1:n.paths){
    NumberOfSOMAproteinsInPathway <- top.ids[4,pathi] %>% as.numeric()              # number of SOMA proteins in reactome
    
    k_proteins_in_class <- NumberOfProteinsSignInResCat
    x_proteins_in_class_n_path <- top.ids[MultestResCati+4,pathi] %>% as.numeric()   # number of proteins in res class in reactome
    m_protein_set <- NumberOfSOMAproteinsInPathway
    n_not_protein_set <- TotalNumberOfSOMAProteins - NumberOfSOMAproteinsInPathway
    
    A <- x_proteins_in_class_n_path
    B <- k_proteins_in_class - x_proteins_in_class_n_path
    C <- m_protein_set - x_proteins_in_class_n_path
    D <- n_not_protein_set - (k_proteins_in_class - x_proteins_in_class_n_path)
    
    # I used this for guidance https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
    
    conTab <- matrix(c(A,B,C,D), 
                     nrow=2, dimnames=list(c("path","npath"), c("sign", "nsign")))
    pvaluestab[pathi,MultestResCati] <- fisher.test(x = conTab)$p.value
    direction_of_fisher_tab[pathi,MultestResCati] <- fisher.test(x = conTab)$estimate
    
  }
}






################## ============== Prepare results for plotting ========== ##############
# Adjust p values within classes
pvaluestab <- apply(pvaluestab,2,p.adjust)

cov.mx <- top.ids.coverage[c(4:7),] %>% as.matrix() %>% t() %>% apply(.,2,as.numeric)
p.cutoff.mx <- ifelse(pvaluestab[,1:4] < 0.05,1,0)

# Multiply so not significant results becomes 0
heatmap.mx <- cov.mx*p.cutoff.mx
# Multiply so significant odds ratios on each side of 1 gets + or -
heatmap.mx <- heatmap.mx*ifelse(direction_of_fisher_tab[,1:4] < 1,0,1)

rownames(heatmap.mx) <- colnames(top.ids.coverage)
colnames(heatmap.mx) <- c("Released to mother","Taken up from mother","Released to fetus","Taken up from fetus")
if(names(table(is.na(heatmap.mx)))){cat("all values in heatmap.mx is na!")}

if(placenta_specific){
  heatmap.mx <- heatmap.mx[,c(1,2)]
  colnames(heatmap.mx) <-  c("Specifically release to mother","Upregulated release to mother")
} else{heatmap.mx <- heatmap.mx[,c(1,2,3,4)]}







######################### ================== Heat map =================== #####################
par(mfrow=c(1,1))
ht <- ComplexHeatmap::Heatmap(heatmap.mx,
                              cluster_rows=FALSE,
                              # row_split = anchestorwords$clusterNames, 
                              # row_title_gp = gpar(fontsize=3),
                              # row_title_rot = 0,
                              row_names_gp = gpar(fontsize = 10),
                              column_names_gp = gpar(fontsize = 10),
                              col = c("white", "yellow", "orange", "red"), 
                              cluster_columns = FALSE,
                              heatmap_legend_param = list(title = "Coverage: % of total Reactome categories",
                                                          direction = "horizontal",
                                                          legend_width = unit(3, "cm"),
                                                          title_position = "leftcenter"))

ComplexHeatmap::draw(ht, 
                     heatmap_legend_side = "bottom", 
                     column_title = paste(), 
                     column_title_gp = gpar(fontsize = 15, fontface = "bold"))





######################### =============== Save heat map to journal =========== ##############
tiff(paste0('SupplementaryFigure',ifelse(placenta_specific,"5_enrheatlay2","4_enrheatlay1"),'.tiff'),
     units="px", width=1500, height=2300, res=300, compression = 'lzw')
ComplexHeatmap::draw(ht, 
                     heatmap_legend_side = "bottom", 
                     column_title = paste(), 
                     column_title_gp = gpar(fontsize = 15, fontface = "bold")); dev.off()