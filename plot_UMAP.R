library(RColorBrewer)
library(ggplot2) 
library(dbscan)

cutoff <- function(num,u_lim,l_lim){
  if(num > u_lim){
    return(u_lim)
  } else if(num < l_lim){
    return(l_lim)
  }else return(num)
}

n_imp = 10
cols <-brewer.pal(3, "Set1")[c(1,2)]
cols2 <-brewer.pal(8, "Set1")[c(4,8)]

original_data <- read.table("pji-blood-extract - cut_c4_fin+CRP_20191108.csv",sep=",",header=T)
original_data <- original_data[rowSums(is.na(original_data)) < 5, ]
original_data_preope <- original_data[original_data$Time_original=="1m_ago",]
original_data_postope <- original_data[original_data$Time_original!="1m_ago",]

all_list <- 1:nrow(original_data_preope)
train_list <- all_list[all_list %% 3 != 0]
test_list <- all_list[all_list %% 3 == 0]

original_data_train <- original_data_preope[train_list,]
original_data_train$dataset <- rep("train",nrow(original_data_train))
original_data_test <- original_data_preope[test_list,]
original_data_test$dataset <- rep("test",nrow(original_data_test))

original_data_postope$dataset <- rep("postope",nrow(original_data_postope))

UMAPdata_train <- read.table("YCU_PJI_preope_train_SRFdist_UMAP_embedding.txt",sep="\t",header=T,row.names=1)
UMAPdata_train <- data.frame(UMAP1=UMAPdata_train[,1],UMAP2=UMAPdata_train[,2],original_data_train)
UMAPdata_test <- read.table("YCU_PJI_preope_test_SRFdist_UMAP_embedding.txt",sep="\t",header=T,row.names=1)
UMAPdata_test <- data.frame(UMAP1=UMAPdata_test[,1],UMAP2=UMAPdata_test[,2],original_data_test)
UMAPdata_postope <- read.table("YCU_PJI_postope_SRFdist_UMAP_embedding.txt",sep="\t",header=T,row.names=1)
UMAPdata_postope <- data.frame(UMAP1=UMAPdata_postope[,1],UMAP2=UMAPdata_postope[,2],original_data_postope)

UMAPdata_preope <- rbind(UMAPdata_train,UMAPdata_test)

UMAPdata_all <- rbind(UMAPdata_preope,UMAPdata_postope)

# res <- dbscan(x = UMAPdata_all[,c("UMAP1","UMAP2")], eps = 0.6, minPts = 4)

res <- kmeans(UMAPdata_all[,c("UMAP1","UMAP2")],2)

clu <- res$cluster
UMAPdata_all$cluster <- factor(clu)

UMAPdata_preope$PJI <- as.factor(UMAPdata_preope$PJI)
UMAPdata_preope$dataset <- factor(UMAPdata_preope$dataset,levels=c("train","test"))
UMAPdata_preope$cluster <- UMAPdata_all$cluster[c(1:nrow(UMAPdata_preope))]

p_norm <- ggplot(data=UMAPdata_all,aes(x=UMAP1,y=UMAP2,color=cluster))+
  geom_point()+
  stat_ellipse(type="norm")

pb <- ggplot_build(p_norm)
ellipse <- pb$data[[2]]
colnames(ellipse) <- c("colour","UMAP1","UMAP2","PANEL","group","size","linetype","alpha")
ellipse1 <- ellipse[ellipse$colour=="#A3A500",]
ellipse2 <- ellipse[ellipse$colour=="#00BF7D",]

p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_preope) +
  geom_point(aes(shape=dataset,fill=PJI),alpha=0.8,size=3,stroke=0.5) +
  scale_shape_manual(values=c(21,23))+
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  #scale_colour_brewer(palette = "Set1") + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  scale_fill_manual(values=c("gray","red"))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)+
  theme(aspect.ratio=1)
ggsave(file = "YCU_PJI_preope_SRFdist_UMAP_diagnosis.pdf", plot = p,width=5,height=4)

continuous_vars <- colnames(original_data)[6:ncol(original_data)]

for (item in continuous_vars){
  UMAPdata_ex <- na.omit(data.frame(UMAPdata_preope[,c("UMAP1","UMAP2","cluster")],UMAPdata_preope[,item,drop=F]))
  colnames(UMAPdata_ex) <- c("UMAP1","UMAP2","cluster","value")

  p <- ggplot(aes(x=cluster, y=value), data=UMAPdata_ex) +
    geom_boxplot(color="black",outlier.shape = NA)+
    geom_jitter(aes(color=cluster),position=position_jitter(0.2),size=2,shape=16)+
    scale_colour_manual(values=cols)+
    theme_classic(base_size = 16)+
    theme(aspect.ratio=2)
    ggsave(file = paste0("YCU_PJI_preope_cluster_boxplot_",item,".pdf"), plot = p,width=4,height=3)


  # u_lim <- quantile(UMAPdata_ex$value,0.95)
  # l_lim <- quantile(UMAPdata_ex$value,0.05)
  # UMAPdata_ex$value <- sapply(UMAPdata_ex$value,u_lim=u_lim,l_lim=l_lim,cutoff)

  if (item == "CRP"){
    UMAPdata_ex$value <- log10(UMAPdata_ex$value)
  }

  p <- ggplot(aes(x=UMAP1, y=UMAP2),data=UMAPdata_ex) +
  geom_point(aes(fill=value),shape=21,alpha=0.8,size=3) +
  scale_fill_distiller(palette = "Spectral") +
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  # scale_colour_manual(values = c("gray","red")) + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  #scale_fill_manual(values=rev(cols[c(1,2)]))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)+
  theme(aspect.ratio=1)
  ggsave(file = paste0("YCU_PJI_preope_SRFdist_UMAP_",item,".pdf"), plot = p,width=5,height=4)
}

antibacteria_bacteria_PJI <- read.table("pji_koukinyaku_plus-bg1_bacteria_PJI.txt",sep="\t",header=T,row.names=1)
antibacteria_bacteria_list <- as.vector(colnames(antibacteria_bacteria_PJI)[c(9:23)])


UMAPdata_PJI <- UMAPdata_preope[UMAPdata_preope$PJI=="1",]
UMAPdata_PJI <- data.frame(UMAPdata_PJI,antibacteria_bacteria_PJI[as.vector(UMAPdata_PJI$ID),])

bac <- c()

for (i in 1:nrow(UMAPdata_PJI)){
  if (!is.na(UMAPdata_PJI$Gram_Positive_cocci[i])){
    bac <- c(bac,"gram+")
  } else if (!is.na(UMAPdata_PJI$Gram_Negative_bacilli[i])){
    bac <- c(bac,"gram-")
  } else {
    bac <- c(bac,"ND")
  }
}

UMAPdata_PJI$bacteria <- factor(bac,levels=c("gram+","gram-","ND"))
UMAPdata_PJI <- UMAPdata_PJI[rev(order(UMAPdata_PJI$bacteria)),]

fisher_mat <- matrix(c(
  nrow(UMAPdata_PJI[UMAPdata_PJI$cluster=="1" & UMAPdata_PJI$bacteria=="gram+",]),
  nrow(UMAPdata_PJI[UMAPdata_PJI$cluster=="1" & UMAPdata_PJI$bacteria!="gram+",]),
  nrow(UMAPdata_PJI[UMAPdata_PJI$cluster=="2" & UMAPdata_PJI$bacteria=="gram+",]),
  nrow(UMAPdata_PJI[UMAPdata_PJI$cluster=="2" & UMAPdata_PJI$bacteria!="gram+",])
  ),ncol=2,byrow=T)

fisher.res <- fisher.test(fisher_mat)

p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_PJI) +
  geom_point(aes(fill=bacteria),shape=21,alpha=0.8,size=3,stroke=0.5) +
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  #scale_colour_manual(values=cols[c(2,3)])+
  # scale_colour_manual(values=cols) +
  scale_fill_manual(values=c(cols2,"gray")) + 
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)+
  theme(aspect.ratio=1)
ggsave(file = paste0("YCU_PJI_SRFdist_UMAP_PJI_bacteria.pdf"), plot = p,width=5,height=4)


bginfo <- read.table("PJI_bginfo.txt",sep="\t",header=T,row.names=1)

UMAPdata_bginfo <- data.frame(UMAPdata_preope,bginfo[as.vector(UMAPdata_preope$ID),])

numerical_vars <- c("height","weight","age","BMI")
categorical_vars <- c("sex","HBP","Diabetes","Smoking","Alcohol")

for (item in numerical_vars){
  UMAPdata_ex <- na.omit(data.frame(UMAPdata_bginfo[,c("UMAP1","UMAP2","cluster")],UMAPdata_bginfo[,item,drop=F]))
  colnames(UMAPdata_ex) <- c("UMAP1","UMAP2","cluster","value")

  p <- ggplot(aes(x=cluster, y=value), data=UMAPdata_ex) +
    geom_boxplot(color="black",outlier.shape = NA)+
    geom_jitter(aes(color=cluster),position=position_jitter(0.2),size=2,shape=16)+
    scale_colour_manual(values=cols)+
    theme_classic(base_size = 16)+
    theme(aspect.ratio=2)
    ggsave(file = paste0("YCU_PJI_preope_cluster_boxplot_",item,".pdf"), plot = p,width=4,height=3)


  # u_lim <- quantile(UMAPdata_ex$value,0.95)
  # l_lim <- quantile(UMAPdata_ex$value,0.05)
  # UMAPdata_ex$value <- sapply(UMAPdata_ex$value,u_lim=u_lim,l_lim=l_lim,cutoff)
  if (item == "CRP"){
    UMAPdata_ex$value <- log10(UMAPdata_ex$value)
  }

  p <- ggplot(aes(x=UMAP1, y=UMAP2),data=UMAPdata_ex) +
  geom_point(aes(color=value),shape=16,alpha=0.8,size=3) +
  scale_colour_distiller(palette = "Spectral") +
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  # scale_colour_manual(values = c("gray","red")) + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  #scale_fill_manual(values=rev(cols[c(1,2)]))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)+
  theme(aspect.ratio=1)
  ggsave(file = paste0("YCU_PJI_preope_SRFdist_UMAP_",item,".pdf"), plot = p,width=4,height=3)
}


for (item in categorical_vars){
  UMAPdata_ex <- na.omit(data.frame(UMAPdata_bginfo[,c("UMAP1","UMAP2","cluster")],UMAPdata_bginfo[,item,drop=F]))
  colnames(UMAPdata_ex) <- c("UMAP1","UMAP2","cluster","value")

  UMAPdata_ex$value <- factor(UMAPdata_ex$value)

  p <- ggplot(aes(x=UMAP1, y=UMAP2),data=UMAPdata_ex) +
  geom_point(aes(fill=value),shape=21,alpha=0.8,size=3) +
  scale_fill_manual(values=c(NA,"black"))+
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  # scale_colour_manual(values = c("gray","red")) + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  #scale_fill_manual(values=rev(cols[c(1,2)]))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)+
  theme(aspect.ratio=1)
  ggsave(file = paste0("YCU_PJI_preope_SRFdist_UMAP_",item,".pdf"), plot = p,width=4,height=3)
}



