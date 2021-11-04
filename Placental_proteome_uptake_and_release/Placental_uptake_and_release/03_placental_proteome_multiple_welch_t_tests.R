###
#
# author: Maren-Helene Høie Degnes
# date: 04.11.21
# purpose: perform welch t-tests to identify the six classes of proteins:
#          - proteins released from placenta to maternal circulation
#          - proteins taken up from maternal circulation to placenta
#          - proteins released from placenta to fetal circulation
#          - proteins taken up from fetal circulation to placenta
#          - proteins upregulated released from placenta to maternal circulation
#          - proteins specifically released from placenta to maternal circulation
#
###




library(haven)
library(tidyverse)
library(scales)
library(openxlsx)





options(stringsAsFactors = FALSE)
rm(list = ls())





###################### ================ Function for paired welch t-tests ============ ###############
PairMulTestD <- function(difftab,difftab2=NA){
  test <- c()
  if(is.na(difftab2)){
    for (i in 1:ncol(difftab)){
      each.test <- t.test(difftab[,i])
      test <- rbind(test, c(each.test$estimate, 
                            each.test$statistic, 
                            each.test$conf.int[1], 
                            each.test$conf.int[2], 
                            each.test$p.value,NA))
    }
    test[,6] <- p.adjust(test[,5], method="BH")
    colnames(test) <- c("mean of differences", 
                        "t statistic", 
                        "lower .95 CI", 
                        "upper .95 CI", 
                        "p value", 
                        "adjusted p")
  } else{  
    
    for (i in 1:ncol(difftab)){
    each.test <- t.test(difftab[,i],difftab2[,i])
    test <- rbind(test, c(each.test$estimate[1],
                          each.test$estimate[2],
                          each.test$statistic, 
                          each.test$conf.int[1], 
                          each.test$conf.int[2], 
                          each.test$p.value,NA))
    }
    colnames(test) <- c("mean of differences1", 
                        "mean of differences1", 
                        "t statistic", 
                        "lower .95 CI", 
                        "upper .95 CI", 
                        "p value", 
                        "adjusted p")}
  test <- test %>% as.data.frame()
  
  
  return(test)
}

  





############## =============== Perform and save results from tests ================ ##############
load("data/lgRFUmedians.RData")

MatNormTest <- PairMulTestD(difftab = adj.m.c.diffs)
FetNormTest <- PairMulTestD(difftab = adj.f.c.diffs)
ArmNormTest <- PairMulTestD(difftab = adj.marm.c.diffs)

secrprot <- which(MatNormTest$`adjusted p`<0.05)

MatPLARMtest <- PairMulTestDD(difftab = adj.m.c.diffs[,secrprot], 
                              difftab2 = adj.marm.c.diffs[,secrprot])
MatNormTest$PLARMpadj <- NA
MatNormTest$PLARMt <- NA
MatNormTest$PLARMpadj[secrprot] <- MatPLARMtest$`adjusted p`
MatNormTest$PLARMt[secrprot] <- MatPLARMtest$`t statistic`

save(MatNormTest,
     FetNormTest,
     ArmNormTest,
     MatPLARMtest,
     file="data/MTT_FBmedians.RData")






############ ====================== p value histogram ============= ############
tiff('supplfig3_phist.tiff',
     units="px", width=2250, height=1200, res=300, compression = 'lzw')


load(paste0("data/MTT_FBmedians.RData"))


b <- 0.025
ylims <- c(0,1000)
par(mfcol=c(1,2),mar=c(4,4,1,0.5), oma=c(1.5,2,2,1))

p <- rep(NA, length(MatNormTest$`p value`))
hist(MatNormTest$`p value`, 
     breaks=seq(0,1,b), 
     main=paste0("A."),
     xlab="",
     ylim = ylims)
# abline(h=qbinom(0.95, length(MatNormTest$`p value`[which(MatNormTest$`t statistic`>=0)]), b), col="orange")
# abline(h=qbinom(1-b*0.05, length(MatNormTest$`p value`[which(MatNormTest$`t statistic`>=0)]), b), col="blue")

hist(FetNormTest$`p value`, 
     breaks=seq(0,1,b), 
     main=paste0("B."),
     xlab="",
     ylim = ylims)
# abline(h=qbinom(0.95, length(MatNormTest$`p value`[which(MatNormTest$`t statistic`<0)]), b), col="orange")
# abline(h=qbinom(1-b*0.05, length(MatNormTest$`p value`[which(MatNormTest$`t statistic`<0)]), b), col="blue")

dev.off()


################# =============== Write excel sheet of all results! ========== ###############
sfx <- "medians"
load(paste0("data/MTT_FB",sfx,".RData"))
load(paste0("data/lgRFU",sfx,".RData"))

padj.limit <- 0.05


MatNormTest$Fullname <- prot$TargetFullName
MatNormTest$Uniprot <- prot$UniProt
MatNormTest$SomaID <- prot$SomaId
MatNormTest$AptName <- prot$AptName


# Prepare secreted to maternal side table
MatNormTest$VADiffInArm <- (((2^(ArmNormTest$`mean of differences`)) - 1) * 100) %>% round(., 2)
MatNormTest$RFUratioarm <- 2^(ArmNormTest$`mean of differences`)
MatNormTest$log2FCarm <- ArmNormTest$`mean of differences`
MatNormTest$adjparm <- ArmNormTest$`adjusted p`
MatNormTest$VADiffInFet <- (((2^(FetNormTest$`mean of differences`)) - 1) * 100) %>% round(., 2)
MatNormTest$RFUratiofet <- 2^(FetNormTest$`mean of differences`)
MatNormTest$log2FCfet <- FetNormTest$`mean of differences`
MatNormTest$adjpfet <- FetNormTest$`adjusted p`


# New column with proteins sign positive difference in mat placenta AND higher than difference in arm!

# Choose proteins secreted to maternal side
spec.idx <- which(MatNormTest$`adjusted p`<0.05
                  & MatNormTest$`t statistic`>0
                  & ((ArmNormTest$`adjusted p`>0.05)
                     | ArmNormTest$`adjusted p`<0.05 & ArmNormTest$`t statistic`<0))

upr.idx <- which(MatNormTest$`adjusted p`<0.05
                 & MatNormTest$`t statistic`>0
                 & MatNormTest$PLARMpadj<0.05
                 & MatNormTest$PLARMt>0)

MatNormTest$MATplupr <- "Not sign plsupr"
MatNormTest$MATplupr[upr.idx] <- "Upregulated release to mother"
MatNormTest$MATplspec <- "Not sign plspec"
MatNormTest$MATplspec[spec.idx] <- "Specific release to mother"

proteinssecretedtomat <- which(MatNormTest$`adjusted p`<0.05
                               & MatNormTest$`t statistic`>0)
proteinstakenupfrommat <- which(MatNormTest$`adjusted p`<0.05
                                & MatNormTest$`t statistic`<0)

MatNormTest$MAT <- "Not secreted or taken up"
MatNormTest$MAT[proteinssecretedtomat] <- "To mother's circulation"
MatNormTest$MAT[proteinstakenupfrommat] <- "From mother's circulation"

MatNormTest$FET <- "Not secreted or taken up"
proteinssecretedtofet <- which(FetNormTest$`adjusted p`<0.05
                               & FetNormTest$`t statistic`>0)
proteinstakenupfromfet <- which(FetNormTest$`adjusted p`<0.05
                                & FetNormTest$`t statistic`<0)

MatNormTest$FET[proteinssecretedtofet] <- "To fetus' circulation"
MatNormTest$FET[proteinstakenupfromfet] <- "From fetus' circulation"

matsecr <- which(MatNormTest$`adjusted p`< padj.limit & MatNormTest$`t statistic`>0)
armsecr <- which(ArmNormTest$`adjusted p`< padj.limit & ArmNormTest$`t statistic`>0)
fetsecr <- which(FetNormTest$`adjusted p`< padj.limit & FetNormTest$`t statistic`>0)
armmatsecr <- intersect(armsecr, matsecr)
matfetsecr <- intersect(fetsecr, matsecr)
fetarmsecr <- intersect(fetsecr, armsecr)

matupt <- which(MatNormTest$`adjusted p`< padj.limit & MatNormTest$`t statistic`<0)
armupt <- which(ArmNormTest$`adjusted p`< padj.limit & ArmNormTest$`t statistic`<0)
fetupt <- which(FetNormTest$`adjusted p`< padj.limit & FetNormTest$`t statistic`<0)
armmatupt <- intersect(matupt, armupt)
matfetupt <- intersect(matupt, fetupt)
fetarmupt <- intersect(armupt, fetupt)

matuptarmsecr <- intersect(armsecr,matupt)
matsecrarmupt <- intersect(matsecr,armupt)
fetuptarmsecr <- intersect(armsecr,fetupt)
fetsecrarmupt <- intersect(fetsecr,armupt)


MatNormTest$secr <- rep("NOT SIGNIF",nrow(MatNormTest))
MatNormTest$secr[matupt] <- "from mat"
MatNormTest$secr[armupt] <- "arm artery protein"
MatNormTest$secr[fetupt] <- "from fet"
MatNormTest$secr[armmatupt] <- "from mat and high artery arm"
MatNormTest$secr[matfetupt] <- "from mat and fet"
MatNormTest$secr[fetarmupt] <- "from fet and high artery arm"
MatNormTest$secr[matsecr] <- "to mat"
MatNormTest$secr[armsecr] <- "arm vein protein"
MatNormTest$secr[fetsecr] <- "to fet"
MatNormTest$secr[armmatsecr] <- "to mat and high vein arm"
MatNormTest$secr[matfetsecr] <- "to mat and fet"
MatNormTest$secr[fetarmsecr] <- "to fet and high vein arm"
MatNormTest$secr[matuptarmsecr] <- "from mat and high vein arm"
MatNormTest$secr[matsecrarmupt] <- "to mat and high artery arm"
MatNormTest$secr[fetuptarmsecr] <- "from fet and high vein arm"
MatNormTest$secr[fetsecrarmupt] <- "to fet and high artery arm"
MatNormTest$secr[intersect(matsecr, fetupt)] <- "from fet to mat"
MatNormTest$secr[intersect(matupt, fetsecr)] <- "from mat to fet"

MatNormTest$secr[intersect(intersect(armupt, matupt), fetupt)] <- "from mat and fet and high artery arm"
MatNormTest$secr[intersect(intersect(matsecr, armupt), fetupt)] <- "from fet to mat and high artery arm"
MatNormTest$secr[intersect(intersect(matupt, armupt), fetsecr)] <- "from mat to fet and high artery arm"
MatNormTest$secr[intersect(intersect(armupt, matsecr), fetsecr)] <- "to mat and fet and high artery arm"
MatNormTest$secr[intersect(intersect(armsecr, matsecr), fetsecr)] <- "to mat and fet and high vein arm"
MatNormTest$secr[intersect(intersect(matupt, armsecr), fetsecr)] <- "from mat to fet and high vein arm"
MatNormTest$secr[intersect(intersect(matsecr, armsecr), fetupt)] <- "from fet to mat and high vein arm"
MatNormTest$secr[intersect(intersect(armsecr, matupt), fetupt)] <- "from mat and fet and high vein arm"



matsec <- MatNormTest
matsec$percent.fold.change <- (((2^(matsec$`mean of differences`)) - 1) * 100) %>% round(.,2)
matsec$RFUratiomat <- 2^(MatNormTest$`mean of differences`)
matsec$log2FCmat <- MatNormTest$`mean of differences`
matsec$CI.95 <- paste0("[", (((2^(matsec$`lower .95 CI`)) - 1) * 100) %>% round(.,2), ",",  (((2^(matsec$`upper .95 CI`)) - 1) * 100) %>% round(.,2), "]")

# mat.in.level <- apply(MatNormPAE, 2, mean)
# matsec$maternalPAE <- (((2^mat.in.level) - 1) * 100)

Mat_ordered <- matsec#[order(matsec$`adjusted p`, decreasing = F),]

table(Mat_ordered$secr)

restab <-   data.frame(Proteinname = Mat_ordered$Fullname, 
                       UniProt = Mat_ordered$Uniprot,
                       SomaID = Mat_ordered$SomaID,
                       AptamerName = Mat_ordered$AptName,
                       "PercFCMat" = Mat_ordered$percent.fold.change,
                       "RFUratioMat" = Mat_ordered$RFUratiomat,
                       "log2FCMat" = Mat_ordered$log2FCmat,
                       "BHAdjPMat" = Mat_ordered$`adjusted p`,
                       "PercFCFet" = Mat_ordered$VADiffInFet,
                       "RFUratioFet" = Mat_ordered$RFUratiofet,
                       "log2FCFet" = Mat_ordered$log2FCfet,
                       "BHAdjPFet" = Mat_ordered$adjpfet,
                       "PercFCArm" = Mat_ordered$VADiffInArm,
                       "RFUratioArm" = Mat_ordered$RFUratioarm,
                       "log2FCArm" = Mat_ordered$log2FCarm,
                       "BHAdjPArm" = Mat_ordered$adjparm,
                       "Transport of protein" = Mat_ordered$secr,
                       "Mothercolumn" = Mat_ordered$MAT,
                       "Fetuscolumn" = Mat_ordered$FET,
                       "Mothersplacentaspec" = Mat_ordered$MATplspec,
                       "Mothersplacentaupr" = Mat_ordered$MATplupr)



write.xlsx(restab,
           file=paste0("SupplementaryTable1", sfx,".xlsx"))

write.xlsx(restab[unique(c(upr.idx,spec.idx)),],
           file=paste0("SupplementaryTable1", sfx,".xlsx"),
           sheet=2)












########################## ===================== Volcano plots ======================= ##########################
xlims <- c(-2,2)
legendcex <- 1.5
numcat <- 10
quantlimperc <- 95


par(mfrow=c(2,3))
plot(ArmNormTest$`mean of differences`,
     -log(ArmNormTest$`p value`,10),
     ylab="-log10(p-value)",
     xlab="Estimated Difference PVE-PAE",
     xlim=xlims)
abline(h=-log10(0.05))
text(1.5,-log10(0.05)+1,"P=0.05")

plot(MatNormTest$`mean of differences`,
     -log(MatNormTest$`p value`,10),
     ylab="-log10(p-value)",
     xlab="Estimated Difference UE-PAE",
     xlim=xlims)
abline(h=-log10(0.05))
text(0.5,-log10(0.05)+1,"P=0.05")

plot(FetNormTest$`mean of differences`,
     -log(FetNormTest$`p value`,10),
     ylab="-log10(p-value)",
     xlab="Estimated Difference VUE-AUE",
     xlim=xlims)
abline(h=-log10(0.05))
text(0.8,-log10(0.05)+1,"P=0.05")




VolcanoPlot(ArmNormTest, 
            themain = "VA-differences maternal arm not background sub", 
            color.vec = GenerateColorvecFromMTT(MatNormTest,
                                                quantilelimperc=quantlimperc,
                                                transparency=0.4,
                                                num.cat=numcat))
VolcanoPlot(MatNormTest, themain = "VA-differences maternal side not background sub", 
            color.vec = GenerateColorvecFromMTT(MatNormTest,
                                                quantilelimperc=quantlimperc,
                                                transparency=0.4, 
                                                num.cat=numcat))
abline(h=-log10(quantile(MatNormTest$`p value`, 
                         (100-quantlimperc)/100)))
text(0.5,-log10(quantile(MatNormTest$`p value`, 
                         (100-quantlimperc)/100))+1,
     paste(quantlimperc, "perc. P=",
           quantile(MatNormTest$`p value`, 
                    (100-quantlimperc)/100) %>% round(6)))
VolcanoPlot(FetNormTest, themain = "VA-differences fetal side not background sub", 
            color.vec = GenerateColorvecFromMTT(FetNormTest,
                                                quantilelimperc=quantlimperc,
                                                transparency=0.4, 
                                                num.cat=numcat))
abline(h=-log10(quantile(FetNormTest$`p value`, 
                         (100-quantlimperc)/100)))
text(0.5,-log10(quantile(FetNormTest$`p value`, 
                         (100-quantlimperc)/100))+1,paste(quantlimperc, "perc. P=", 
                                                          quantile(FetNormTest$`p value`, 
                                                                   (100-quantlimperc)/100) %>% round(6)))


