library(lawstat)
library(corrplot)
library(readr)
library(ggplot2)

df <- read_csv("pji-blood-extract - cut_c4_fin+CRP_20191108.csv")
df <- df[df$Time==1,]

continuous_vars <- colnames(df[,6:ncol(df)])

continuous_res <- c()

for (v in continuous_vars){
	df_ind <- na.omit(data.frame(PJI=df$PJI,var=df[,v]))
	colnames(df_ind) <- c("PJI","var")

	res <- brunner.munzel.test(df_ind[df_ind$PJI==0,"var"],df_ind[df_ind$PJI==1,"var"])
	res2 <- t.test(x=df_ind[df_ind$PJI==0,"var"],y=df_ind[df_ind$PJI==1,"var"],var.equal=F,paired=F)
	continuous_res <- rbind(continuous_res,c(v,res$statistic,res$p.value,res2$statistic,res2$p.value))
}

continuous_res <- data.frame(continuous_res)
colnames(continuous_res) <- c("variable","BM.statistic","BM.p.value","Welch.statistic","Welch.p.value")
write.table(continuous_res,file="PJI_stattest_continuous_vars.txt",sep="\t",quote=F,row.names=F)

df_ex <- na.omit(df[,c("PJI","CRP","hemoglobin")])

df1 <- df_ex[df_ex$CRP < 0.3 & df_ex$hemoglobin < 11,]
df2 <- df_ex[df_ex$CRP < 0.3 & df_ex$hemoglobin >= 11,]
df3 <- df_ex[df_ex$CRP >= 0.3 & df_ex$hemoglobin < 11,]
df4 <- df_ex[df_ex$CRP >= 0.3 & df_ex$hemoglobin >= 11,]

df$PJI <- as.factor(df$PJI)

ggplot(data=df,aes(x=log10(CRP),y=hemoglobin)) +
	geom_hline(yintercept = 11, linetype = "dashed") +
	geom_vline(xintercept = log10(0.3), linetype = "dashed") +
	geom_point(aes(fill=PJI),shape=21,alpha=0.8,size=3,stroke=0.5) +
	scale_fill_manual(values=c("gray","red"))+
	ylim(9,17)+
	theme_classic(base_size = 16)+
 	theme(aspect.ratio=1)
 ggsave(file = "YCU_PJI_preope_CRP_Hb.pdf",width=5,height=4)

