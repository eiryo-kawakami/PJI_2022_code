library(pROC)
library(RColorBrewer)
library(matrixStats)


mycols <-brewer.pal(5, "Set1")

prob_group1 <- read.table("YCU_PJI_preope_RFprob_calibrated_isotonic.txt",sep="\t",header=T,row.names=1)

roc.res.group1 <- roc(prob_group1$PJI, rowMedians(as.matrix(prob_group1[,c(4:ncol(prob_group1))])))

print(pROC::auc(roc.res.group1))

rocplot <- "YCU_PJI_RF_rocplot.pdf"
pdf(rocplot,useDingbats=FALSE)
plot.roc(roc.res.group1,col="black")
dev.off()

prob_group1 <- read.table("YCU_PJI_preope_LRprob.txt",sep="\t",header=T,row.names=1)

roc.res.group1 <- roc(prob_group1$PJI, rowMedians(as.matrix(prob_group1[,c(4:ncol(prob_group1))])))

print(pROC::auc(roc.res.group1))

rocplot <- "YCU_PJI_LR_rocplot.pdf"
pdf(rocplot,useDingbats=FALSE)
plot.roc(roc.res.group1,col="black")
dev.off()

prob_group1 <- read.table("YCU_PJI_preope_XGBprob.txt",sep="\t",header=T,row.names=1)

roc.res.group1 <- roc(prob_group1$PJI, rowMedians(as.matrix(prob_group1[,c(4:ncol(prob_group1))])))

print(pROC::auc(roc.res.group1))

rocplot <- "YCU_PJI_XGB_rocplot.pdf"
pdf(rocplot,useDingbats=FALSE)
plot.roc(roc.res.group1,col="black")
dev.off()
