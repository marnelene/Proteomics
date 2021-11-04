###
#
# author: Maren-Helene Høie Degnes
# date: 04.11.21
# purpose: Investigate the gestational patterns of upregulated and/or specifically released proteins from placenta to mother 
#
###

library(tidyverse)
library(Hotelling)
library(lme4)
library(splines)
library(dendextend)
library(ape)
library(RColorBrewer)
library(cluster)
library(scales)
library(openxlsx)



###Some decisions
plot_VAarmpl <- F #readline(prompt = "Do you want to plot the differences between V-A in placenta against arm for all uprand/orspec proteins? T/F   ")
mod_onlyselected <- F #readline(prompt = "Do you want to model only selected 6 proteins? T/F  ")
plotpat_eachmodelling <- T #readline(prompt = "Do you want to plot each protein for each modelling? T/F  ")
mod_only_changing <- F





############### ============= load dataset and subset upregulated and/or specific proteins ============== ##################
sfx <- "medians"
load(paste0("data/lgRFU",sfx,".RData"))
load(paste0("MTT_FB",sfx,".RData"))
samp.ordering <- order(samp$GAWeeks);    samprot <- samprot[samp.ordering,];    samp <- samp[samp.ordering,]  # order by gaweeks for plotting later

# upregulated proteins are 
#             1. higher secreted from placenta to mother than in arm
#             
upr.idx <- which(MatNormTest$`adjusted p`<0.05
                 & MatNormTest$`t statistic`>0
                 & MatNormTest$PLARMpadj<0.05
                 & MatNormTest$PLARMt>0)

# specific proteins are 
#             1. secreted from placenta to mother
#             2. AND difference in arm should be 0 or negative
spec.idx <- which(MatNormTest$`adjusted p`<0.05
                  & MatNormTest$`t statistic`>0
                  & ((ArmNormTest$`p value`>0.05)
                     | ArmNormTest$`t statistic`<0))

# Save the three different categories 
spec_but_not_upr <- setdiff(spec.idx,upr.idx)
spec_and_upr <- intersect(spec.idx,upr.idx)
upr_but_not_spec <- setdiff(upr.idx,spec.idx)
all.pr <- c(spec_but_not_upr,spec_and_upr,upr_but_not_spec) %>% unique()

# Save information about selected proteins for short access later
prodclass.gen <- ifelse(1:nrow(MatNormTest)%in%spec_but_not_upr,1,ifelse(1:nrow(MatNormTest)%in%upr_but_not_spec,3,2))
prodclass <- prodclass.gen[all.pr]
prodclass.protnames <- prot$TargetFullName[all.pr]

if(plot_VAarmpl){## Plot mat diff and arm diff within categories
    par(mfrow=c(1,1))
  plot(c(),ylim=c(-0.3,0.8),xlim=c(0,101),ylab="log2RFU V-A",xlab="protein number");  abline(h=0)
  points(MatNormTest$`mean of differences`[all.pr],col=prodclass+1,pch=16)
  points(ArmNormTest$`mean of differences`[all.pr],col=alpha(prodclass+1,0.4),pch=17)
  legend("topright",legend=c("only specific","specific and upr","only upr","placenta V-A","arm V-A"),col=c(2:4,1,1),pch=c(15,15,15,16,17),bty="n")
}





############## ================== Write excel sheet with information about these proteins ============ #########

uprspec.df <- data.frame(Proteinname=prot$TargetFullName[all.pr],
                         Uniprot=prot$UniProt[all.pr],
                         SOMAID=prot$SomaId[all.pr],
                         AptamerName=prot$AptName[all.pr],
                         "% RFU Mat"= (((2^(MatNormTest$`mean of differences`[all.pr])) - 1) * 100) %>% round(., 2),
                         "log2difference Mat"=MatNormTest$`mean of differences`[all.pr],
                         "FDRcorrP Mat"=MatNormTest$`adjusted p`[all.pr],
                         "% RFU Arm"= (((2^(ArmNormTest$`mean of differences`[all.pr])) - 1) * 100) %>% round(., 2),
                         "log2difference Arm"=ArmNormTest$`mean of differences`[all.pr],
                         "FDRcorrP Arm"=ArmNormTest$`adjusted p`[all.pr],
                         "FDRcorrP Mat-Arm"=MatNormTest$PLARMpadj[all.pr],
                         "Upr and/or spec"=ifelse(prodclass==1,"only spec",ifelse(prodclass==2,"spec and upr","upr")))

write.xlsx(uprspec.df,
           file="SupplementaryTable1",sheet=2)

mean_of_log2RFU_differences <- MatNormTest$`mean of differences`

par(mfrow=c(1,3))
hist(mean_of_log2RFU_differences,breaks=100)
hist(((2^mean_of_log2RFU_differences)),breaks=100)
hist(((2^mean_of_log2RFU_differences-1)*100),breaks=100)




################### ======================= Estimate and plot development curve ========== ##########
# Create containers needed inside for-loop
ST.pred.mx <- c()
pred.mx <- c()
models <- list()
modelsi <- 1
pdfpfx <- "UprSpec_dev"
prot.indeces <- all.pr
bsreslist <- list()   # container for testing coefficients in mm

# Save objects 
normal_longit.idx <- which(!is.na(samp$TimePoint) & samp$Control==1)
normal_4VPVE.idx <- which(samp$SampleGroup=="PVE" & samp$Control==1)
ST.ID <- samp$ID[normal_longit.idx]                          # ID for only stork cohort
meanTP <- tapply(samp$GAWeeks[normal_longit.idx] %>% as.numeric(),samp$TimePoint[normal_longit.idx],mean)


if(mod_onlyselected){  # Shoudl only selected be modelled and plotted? See top at the script 
  prot.indeces <- all.pr[c(which(prot$TargetFullName[all.pr]=="Placenta growth factor"),
                           which(prot$TargetFullName[all.pr]=="Pleiotrophin"),
                           which(prot$TargetFullName[all.pr]=="Secreted frizzled-related protein 3" ),
                           which(prot$TargetFullName[all.pr]=="Tissue factor pathway inhibitor 2" ),
                           58,
                           which(prot$TargetFullName[all.pr]=="Chorionic somatomammotropin hormone"))]
  
  tiff('Suppl_6dev.tiff',
       units="px", width=2250, height=1200, res=300, compression = 'lzw',family = "Arial")
  layout(mat=matrix(c(1,1,2,2,3,3,4,5,5,6,6,7,7,8),nrow=2,byrow=T))
  # par(mar=rep(4,4))
  
} else{pdf('SupplementaryFigure7_UrSpecdevs.pdf',
           width=12.50, height=7.00);   par(mfrow=c(2,3))}


each_fifth_plot_of_six_plots <- seq(5,length(prot.indeces),by=6)
each_first_plot_of_six_plots <- seq(1,length(prot.indeces),by=6)

# Model and plot all proteins indexed inside prot.indeces
for (prnum in 1:length(prot.indeces)){
  par(mar=c(3,3,2,0.5),xpd=T,cex.main=0.9,ps=8)
  pr <- prot.indeces[prnum]
  # set data
  proteinmeas <- samprot[c(normal_longit.idx,normal_4VPVE.idx),pr]
  ST.proteinmeas <- samprot[c(normal_longit.idx),pr]
  proteinname <- prot$TargetFullName[pr]
  uniprot <- prot$UniProt[pr]
  aptname <- prot$AptName[pr]
  
  # patient plots
  patcol <- alpha("darkgrey",0.4)
  patidx <- lapply(unique(ST.ID),function(x){which(ST.ID%in%x)})
  if(plotpat_eachmodelling){
    GAWeeks <- c(samp$GAWeeks[normal_longit.idx],samp$GAWeeks[normal_4VPVE.idx]) %>% as.numeric()
    plot(GAWeeks[patidx[[1]]], proteinmeas[patidx[[1]]], 
         type="l",
         ylab="",xlab="",
         ylim=c(min(proteinmeas), max(proteinmeas)), 
         xlim=c(min(GAWeeks),max(GAWeeks)),
         main=paste(proteinname),
         lwd=0.7)
    mtext(side=3, line=0, cex=0.5, paste("uniprot:",uniprot, "aptname:",aptname))
    for (pat in 2:length(ST.ID %>% unique())){points(GAWeeks[patidx[[pat]]],
                                                     proteinmeas[patidx[[pat]]],
                                                     type="l",
                                                     lwd=0.7)}
    Timepoints <- c(samp$TimePoint[normal_longit.idx],rep(4,length(normal_4VPVE.idx)))
    points(GAWeeks[which(Timepoints==4)],
           proteinmeas[which(Timepoints==4)],
           pch=3,
           lwd=0.7,
           cex=0.7)
  }
  

  ####  mixed modelling for stork + 4V samples
  x <- samp$GAWeeks[c(normal_longit.idx,normal_4VPVE.idx)] %>% as.numeric()
  y <- proteinmeas
  two.cohorts.meanTP <- c(meanTP, mean(samp$GAWeeks[normal_4VPVE.idx]))
  two.cohorts.ID <- samp$ID[c(normal_longit.idx,normal_4VPVE.idx)]
  bsplines <- bs(x, degree=3,knots=two.cohorts.meanTP[c(2)], Boundary.knots =two.cohorts.meanTP[c(1,4)])
  bsplines.df <- bsplines %>% as.data.frame()
  colnames(bsplines.df) <- paste0("bs",1:4)
  mmfit <- lmer(y ~ 1 + bs1 + bs2 + bs3 + bs4 + (1| two.cohorts.ID),
                data=cbind(bsplines.df,two.cohorts.ID,y),REML=F)            #x = ST.GAWeeks, y = observed proteinmeas
  bsreslist[[prnum]] <- lmerTest::as_lmerModLmerTest(mmfit) %>% summary()
  xnew <- seq(two.cohorts.meanTP[1],two.cohorts.meanTP[4],length.out = 100)
  xnew.df <- predict(bsplines, xnew) %>% as.data.frame()
  colnames(xnew.df) <- paste0("bs",1:4)
  newdata <- cbind(xnew.df)
  ynew <- predict(mmfit, newdata=newdata, re.form=~0)                     #ynew = predicted proteinmeas at xnew
  lines(xnew, ynew, type="l",col=alpha(4,0.4),lwd=2.5,lty=2)                      #add a line going through the new predicted protein expression measurements based on the B-splines
  ST.pred.mx <- rbind(ST.pred.mx,ynew)
  
  
  ####  mixed modelling for only stork samples
  x <- samp$GAWeeks[normal_longit.idx] %>% as.numeric()
  y <- ST.proteinmeas
  bsplines <- bs(x, degree=3,knots=meanTP[c(2)], Boundary.knots =meanTP[c(1,3)])
  bsplines.df <- bsplines %>% as.data.frame()
  colnames(bsplines.df) <- paste0("bs",1:4)
  mmfit <- lmer(y ~ 1 + bs1 + bs2 + bs3 + bs4 + (1| ST.ID),
                data=cbind(bsplines.df,ST.ID,y),REML=F)            #x = ST.GAWeeks, y = observed proteinmeas
  bsreslist[[prnum]] <- lmerTest::as_lmerModLmerTest(mmfit) %>% summary()
  xnew <- seq(meanTP[1],meanTP[3],length.out = 100)
  xnew.df <- predict(bsplines, xnew) %>% as.data.frame()
  colnames(xnew.df) <- paste0("bs",1:4)
  newdata <- cbind(xnew.df)
  ynew <- predict(mmfit, newdata=newdata, re.form=~0)                     #ynew = predicted proteinmeas at xnew
  lines(xnew, ynew, type="l",col=3,lwd=2.5)                      #add a line going through the new predicted protein expression measurements based on the B-splines
  ST.pred.mx <- rbind(ST.pred.mx,ynew)
  

  start.point <- mmfit@beta[1]

  #Set legend
  if(mod_onlyselected){
    if(pr==3929){mtext("--Gestational weeks--",1,line=2,cex=0.7)}
    if(pr==2505){mtext("--log2 RFU--",2,line=2,at=8,cex=0.7)}

    if(prnum==3){
      plot(NA ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0,1), ylim=0:1)

      legend(-0.9,1,legend=c("STORK",
                             "4 vessel",
                             "Mean curve STORK",
                             "Mean curve STORK + 4V"),
             lwd=c(1,NA,1.1,1.1),
             pch=c(NA,3,NA,NA,NA),
             col=c("black","black",3,alpha(4,0.4)),
             bty="n",
             lty=c(1,1,1,2))
    }
  }

  else {      legend("topleft",legend=c("STORK","4 vessel",
                                       "Estimated mean curve STORK",
                                       "Estimated mean curve STORK + 4V"),
                       lwd=c(0,NA,1,1),pch=c(NA,3,NA,NA),
                       col=c("black","black",3,alpha(4,0.4)),bty="n",cex=1,
                       lty=c(1,1,1,2))
    if(prnum %in% each_fifth_plot_of_six_plots){mtext("--Gestational weeks--",1,line=2,cex=0.7)}
    if(prnum %in% each_first_plot_of_six_plots){mtext("--log2 RFU--",2,line=2,at=start.point,cex=0.7)}
    
    }
}; dev.off()







################ =============== Investigate only proteins that actually change in all splines ============== #########

if(mod_only_changing){
  spl_est_pvalues <- sapply(bsreslist, function(x){x$coefficients[,5]})  # get p-values
  spl_est_adjpvalues <- apply(spl_est_pvalues, 2, p.adjust)             # adjust p-values by BH, each spline estimate wise
  
  changing_idxvec_idx <- c()           # save indeces for proteins coming out as not changing at one of splines estimates
  nonch.idx <- c()                     # save indeces for proteins coming out as changing at all of splines estimates
  for ( i in 1:ncol(spl_est_adjpvalues)){
    if(spl_est_adjpvalues[2,i]  < 0.05 |    
       spl_est_adjpvalues[3,i]  < 0.05 |  
       spl_est_adjpvalues[4,i]  < 0.05 | 
       spl_est_adjpvalues[5,i]  < 0.05 ){ changing_idxvec_idx <- c(changing_idxvec_idx,i)
    } else{nonch.idx <- c(nonch.idx,i)}
  }
  
  # plot non changing proteins to see if this analysis makes sense
  par(mfrow=c(2,3))
  for (prnum in 1:length(nonch.idx)){
    par(mar=c(3,3,2,0.5),xpd=T,cex.main=0.9,ps=8)
    pr <- prot.indeces[nonch.idx][prnum]
    # set data
    proteinmeas <- samprot[c(normal_longit.idx,normal_4VPVE.idx),pr]
    ST.proteinmeas <- samprot[c(normal_longit.idx),pr]
    proteinname <- prot$TargetFullName[pr]
    uniprot <- prot$UniProt[pr]
    aptname <- prot$AptName[pr]
    
    # patient plots
    patcol <- alpha("darkgrey",0.4)
    patidx <- lapply(unique(ST.ID),function(x){which(ST.ID%in%x)})
    if(T){
      GAWeeks <- c(samp$GAWeeks[normal_longit.idx],samp$GAWeeks[normal_4VPVE.idx]) %>% as.numeric()
      plot(GAWeeks[patidx[[1]]], proteinmeas[patidx[[1]]], 
           type="l",
           ylab="",xlab="",
           ylim=c(min(proteinmeas), max(proteinmeas)), 
           xlim=c(min(GAWeeks),max(GAWeeks)),
           main=paste(proteinname),
           lwd=0.7)
      mtext(side=3, line=0, cex=0.5, paste("uniprot:",uniprot, "aptname:",aptname))
      for (pat in 2:length(ST.ID %>% unique())){points(GAWeeks[patidx[[pat]]],proteinmeas[patidx[[pat]]],type="l",lwd=0.7)}
      Timepoints <- c(samp$TimePoint[normal_longit.idx],rep(4,length(normal_4VPVE.idx)))
      points(GAWeeks[which(Timepoints==4)],proteinmeas[which(Timepoints==4)],pch=3,lwd=0.7,cex=0.7)
      xnew <- seq(meanTP[1],meanTP[3],length.out = 100)
      ynew <- ST.pred.mx[nonch.idx[prnum],]
      lines(xnew,ynew,col=4,lwd=2)
    }
  }
  # sub set to focus only on changing proteins
  ST.pred.mx <- ST.pred.mx[changing_idxvec_idx,]
  prot.indeces <- prot.indeces[changing_idxvec_idx]
  prodclass <- prodclass.gen[all.pr[changing_idxvec_idx]]
  prodclass.protnames <- prot$TargetFullName[all.pr[changing_idxvec_idx]]
}








############## ================= Cluster predicted mean curves ==================== #################

# Create dendrogram objects
cs <- 1-(ST.pred.mx %>% t() %>% cor(.,method="spearman"))#; xnew <- seq(min(GAWeeks),max(GAWeeks),length.out = 100) ; lastt <- 4
csd <- cs %>% as.dist(); hc <- csd %>% hclust(); hc$labels <- prot$TargetFullName[prot.indeces]
dend <- hc %>% as.dendrogram()

# Detect optimal number of clusters above two clusters
silh.means <- c()
for(ncl in 3:10){silh.means <- c(silh.means,summary(silhouette(cutree(hc,k=ncl),csd))$avg.width)}
par(mfrow=c(1,1),mar=(rep(4,4)));plot(3:10,silh.means,type="b",
                                      main="Mean silhouette widths for each number of clusters 3:10",
                                      ylab="Mean silhouette widths")
optimal.num.cl <- grep(max(silh.means),silh.means) + 2

# Cut tree by a cutoff and save colors for plotting
clusMember <- cutree(dend,optimal.num.cl)
allcolors <- c(brewer.pal(n=8,name="Dark2"),"darkcyan","firebrick1")
clusColors <- allcolors[1:optimal.num.cl]
labels_colors(dend) <- clusColors[clusMember][order.dendrogram(dend)]








################==================== Combined plot of dendrogram and scaled developments in tiff ==================###############

tiff('Figure4_dendDev.tiff',
     units="px", width=2250, height=1700, res=300, compression = 'lzw',family = "Arial")





# Set all sizes and the plot window layout
dotcex <- 0.5
leavescex <- 1
axisvaluesdendcex <- 1
axislabelscex <- 1
plotmaincex <- 1
legendcex <- 1
linethickness <- 1
xlabelcex <- 1
layout(mat=matrix(c(rep(1,optimal.num.cl*2),2:(optimal.num.cl+1)),nrow=3,byrow=TRUE))
par(mar=c(16.5,3,0.4,0),
    ps=8,
    cex.axis=axisvaluesdendcex,
    cex.lab=axislabelscex,
    cex.main=plotmaincex,
    mgp=c(3,0.5,0.5))

# connect protein categories (upr/spec) to the dendogram object 
idx.vec <- c()
for (protname in labels(dend)){idx.vec <- c(idx.vec,which(prodclass.protnames==protname))}

# Adjust dendrogram settings
dend <- dend %>% set("leaves_pch", 16) %>%  # node point type
  set("leaves_cex", dotcex) %>%  # node point size
  set("leaves_col", prodclass[idx.vec]) %>% # node point color
  set("labels_cex",1) %>%
  set("branches_lwd", c(0.5,0.5,0.5))

# Create cluster name vector
clusternames <- paste("Cluster", c(2,6,5,1,4,3) #c(6, 1, 3, 4, 2,5)
                      )
# plot dendrogram
plot(dend,edgePar = list(t.col=clusColors),
     ylab="",yaxt="n",lwd=0.2#,ylim=c(0,0.001)
     )
mtext("Spearman correlation",side=2,line=2,outer=F,cex=xlabelcex,at=1)
legend(0,2,col=c(1,2,3),pch=16,
       legend=c("Specific","Specific and Upregulated","Upregulated"),
       bty="n",cex=legendcex)
axis(2,at=seq(0,2,by=0.5),labels=seq(1,-1,by=-0.5))
abline(h=0.35,lty=2,lwd=0.5)

# Settings for development curves beneath dendrogram
par(mar=c(1,1,1,1))
clusterdata <- ST.pred.mx
x <- seq(meanTP[1],meanTP[3],length.out = 100)



firstplot <- T
for (cluster in c(4,1,6,5,3,2) #c(2,5,3,4,6,1)
     ){
  clustercol <- clusColors[cluster]
  clus <- which(clusMember==cluster)
  
  clusterdata[clus,] <- clusterdata[clus,] + (mean(clusterdata[clus,1]) - clusterdata[clus,1])
  clusterdata[clus,] <- clusterdata[clus,] - mean(clusterdata[clus,1])
  
  if(!firstplot){
    if(cluster==4){mtext("--Gestational age in weeks--",side=1,line=2,outer=F,at=33,cex=xlabelcex)}
    par(mar=c(3,2,1,0))
    plot(x,clusterdata[clus[1],],type="l",
         ylim=c(min(clusterdata[clus,]),max(clusterdata[clus,])),col=clustercol,
         xlab="",ylab="",
         lwd=linethickness,
         main=clusternames[cluster])
  }else{
    par(mar=c(3,3,1,0))
    plot(x,clusterdata[clus[1],],type="l",
         ylim=c(min(clusterdata[clus,]),max(clusterdata[clus,])),col=clustercol,
         xlab="",ylab="",lwd=linethickness,
         main=clusternames[cluster])
    mtext("Scaled log2 RFU",side=2,line=2,outer=F,cex=xlabelcex)
  }
  # abline(h=16,col="grey");abline(h=14,col="grey");abline(h=12,col="grey");abline(h=10,col="grey")
  if (length(clus)>1){
    for (protein in clus[2:length(clus)]){
      lines(x,clusterdata[protein,],col=clustercol,lwd=linethickness)
    }
  }
  firstplot <- F
}




dev.off()





