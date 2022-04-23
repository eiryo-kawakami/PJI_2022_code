library(ggplot2)
library(matrixStats)
library(tidyr)


RFprob <- read.table("YCU_PJI_preope_RFprob_calibrated_isotonic.txt",sep="\t",header=TRUE)
RFprob$prob <- rowMedians(as.matrix(RFprob[,c(5:ncol(RFprob))]))
RFprob$PJI <- factor(RFprob$PJI)

LRprob <- read.table("YCU_PJI_preope_LRprob.txt",sep="\t",header=TRUE)
LRprob$prob <- rowMedians(as.matrix(LRprob[,c(5:ncol(LRprob))]))
LRprob$PJI <- factor(LRprob$PJI)

XGBprob <- read.table("YCU_PJI_preope_XGBprob.txt",sep="\t",header=TRUE)
XGBprob$prob <- rowMedians(as.matrix(XGBprob[,c(5:ncol(XGBprob))]))
XGBprob$PJI <- factor(XGBprob$PJI)

prob_summary1 <- data.frame(ID=RFprob$ID,PJI=RFprob$PJI,RF=RFprob$prob,LR=LRprob$prob)
prob_summary2 <- tidyr::gather(prob_summary1,key=method,value=prob,"RF","LR")

ggplot(data=prob_summary2,aes(x=method,y=prob,fill=PJI)) +
	geom_line(aes(color=PJI,group=ID),linetype=2)+
	geom_dotplot(binaxis = "y",stackdir="center",binwidth=0.02,dotsize=1) +
	scale_fill_manual(values=c("gray","red"))+
	scale_color_manual(values=c("gray","red"))+
	ylim(0,1)+
	theme_classic()
ggsave("YCU_PJI_preope_prob.pdf",width=3,height=4,useDingbats=FALSE)

ggplot(data=prob_summary1,aes(x=LR,y=RF,fill=PJI)) +
	geom_abline(slope=1,intercept=0,linetype=2)+
	geom_point(shape=21) +
	scale_fill_manual(values=c("gray","red"))+
	xlim(0,1)+
	ylim(0,1)+
	theme_classic()+
	theme(aspect.ratio=1)
ggsave("YCU_PJI_preope_prob_LRvsRF.pdf",width=5,height=4,useDingbats=FALSE)

