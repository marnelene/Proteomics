###
#
# author: Maren-Helene Høie Degnes
# date: 04.11.21
# purpose: adjust differences between pairs of vessels for water shift and possible other biases using median normalization
#
###

library(tidyverse)
library(haven)
library(pls)




rm(list=ls())
load(file="data/lgRFU.RData")





############## ================ Choose patients from 4V and create new variables ============ ############
UE.samp <- samp[which(samp$SampleGroup=="UE"),]
UE.samprot <- samprot[which(samp$SampleGroup=="UE"),]

PAE.samp <- samp[which(samp$SampleGroup=="PAE"),]
PAE.samprot <- samprot[which(samp$SampleGroup=="PAE"),]

PVE.samp <- samp[which(samp$SampleGroup=="PVE"),]
PVE.samprot <- samprot[which(samp$SampleGroup=="PVE"),]

AUE.samp <- samp[which(samp$SampleGroup=="AUE"),]
AUE.samprot <- samprot[which(samp$SampleGroup=="AUE"),]

VUE.samp <- samp[which(samp$SampleGroup=="VUE"),]
VUE.samprot <- samprot[which(samp$SampleGroup=="VUE"),]






######### ========= Check that patients have the same order so they can be compared and normalized pairwise ======== #######
(UE.samp$ID == PAE.samp$ID) %>% table()
(VUE.samp$ID == AUE.samp$ID) %>% table()
(PVE.samp$ID == PAE.samp$ID) %>% table() # ohno -> # PVE and PAE includes different patients!!
# Find overlap patient ids between PVE and PAE
PVEPAE.patients <- intersect(PVE.samp$ID, PAE.samp$ID)
(PVE.samp$ID[which(PVE.samp$ID  %in% PVEPAE.patients)] == PAE.samp$ID[which(PAE.samp$ID  %in% PVEPAE.patients)]) %>% table()






############ =============== Calculate differences between pairs of vessels ============= ###########
m.c.diffs <- UE.samprot-PAE.samprot
adj.m.c.diffs <- m.c.diffs
for (patient in 1:nrow(m.c.diffs)){ 
  adj.m.c.diffs[patient,] <- m.c.diffs[patient,] - median(m.c.diffs[patient,]) 
}

f.c.diffs <- VUE.samprot-AUE.samprot
adj.f.c.diffs <- f.c.diffs
for (patient in 1:nrow(f.c.diffs)){ 
  adj.f.c.diffs[patient,] <- f.c.diffs[patient,] - median(f.c.diffs[patient,]) 
}


marm.c.diffs <- PVE.samprot[which(PVE.samp$ID  %in% PVEPAE.patients),] - PAE.samprot[which(PAE.samp$ID  %in% PVEPAE.patients),]
adj.marm.c.diffs <- marm.c.diffs
for (patient in 1:nrow(marm.c.diffs)){ 
  adj.marm.c.diffs[patient,] <- marm.c.diffs[patient,] - median(marm.c.diffs[patient,]) 
}






############# ============== control! ============== ################
(prot$AptName == colnames(samprot)) %>% table()
(prot$AptName == colnames(UE.samprot)) %>% table()
(prot$AptName == colnames(adj.m.c.diffs)) %>% table()
(prot$AptName == colnames(m.c.diffs)) %>% table()
(prot$AptName == colnames(adj.marm.c.diffs)) %>% table()
(prot$AptName == colnames(marm.c.diffs)) %>% table()
(prot$AptName == colnames(adj.f.c.diffs)) %>% table()
(prot$AptName == colnames(f.c.diffs)) %>% table()
(rownames(samp) == rownames(samprot)) %>% table()
(apply(adj.m.c.diffs,1,median) == 0) %>% table()
(apply(adj.marm.c.diffs,1,median) == 0) %>% table()
(apply(adj.f.c.diffs,1,median) == 0) %>% table()






############## ============== Write r data file of datasets ============ ############
save(samprot,
     samp,
     prot,
     adj.m.c.diffs,
     adj.marm.c.diffs,
     adj.f.c.diffs,
     file="data/lgRFUmedians.RData")
