library(RColorBrewer)
library(ggplot2) 
library(dbscan)

cutoff <- function(num,u_lim){
  if(num > u_lim){
    return(u_lim)
  }else return(num)
}

n_imp = 10
cols <-brewer.pal(3, "Set1")[c(2,1)]

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

# res <- dbscan(x = UMAPdata_all[,c("UMAP1","UMAP2")], eps = 0.8, minPts = 10)

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

UMAPdata_all$PJI <- as.factor(UMAPdata_all$PJI)
UMAPdata_all$Time_original <- factor(UMAPdata_all$Time_original, levels=c("1m_ago","3m_later","6m_later","12m_later"))

antibacteria_bacteria_PJI <- read.table("pji_koukinyaku_plus-bg1_bacteria_PJI.txt",sep="\t",header=T,row.names=1)

UMAPdata_all <- data.frame(UMAPdata_all,antibacteria_bacteria_PJI[as.vector(UMAPdata_all$ID),])

PJI_cluster <- c()

for (i in 1:nrow(UMAPdata_all)){
  if (nrow(UMAPdata_preope[UMAPdata_preope$ID==UMAPdata_all$ID[i],]) > 0){
    if (UMAPdata_all$PJI[i]=="1"){
      if (UMAPdata_preope[UMAPdata_preope$ID==UMAPdata_all$ID[i],"cluster"][[1]]=="1"){
        PJI_cluster <- c(PJI_cluster,"PJI_cluster1")
      } else {
        PJI_cluster <- c(PJI_cluster,"PJI_cluster2")
      }
    } else {
      if (UMAPdata_preope[UMAPdata_preope$ID==UMAPdata_all$ID[i],"cluster"][[1]]=="1"){
        PJI_cluster <- c(PJI_cluster,"nonPJI_cluster1")
      } else {
        PJI_cluster <- c(PJI_cluster,"nonPJI_cluster2")
      }
    }
  } else {
    PJI_cluster <- c(PJI_cluster,NA)
  }
}

# PJI_cluster <- c()

# for (i in 1:nrow(UMAPdata_all)){
#   if (nrow(UMAPdata_preope[UMAPdata_preope$ID==UMAPdata_all$ID[i],]) > 0){
#     if (UMAPdata_all$PJI[i]=="1"){
#       if (!is.na(UMAPdata_all$Gram_Positive_cocci[i])) {
#         PJI_cluster <- c(PJI_cluster,"PJI_gram+")
#       } else {
#         PJI_cluster <- c(PJI_cluster,"PJI_wo_gram+")
#       }
#     } else {
#       PJI_cluster <- c(PJI_cluster,"nonPJI")
#     }
#   } else {
#     PJI_cluster <- c(PJI_cluster,NA)
#   }
# }

# UMAPdata_all$PJI_cluster <- factor(PJI_cluster,levels=c("PJI_gram+","PJI_wo_gram+","nonPJI"))

UMAPdata_all$PJI_cluster <- factor(PJI_cluster,levels=c("PJI_cluster1","PJI_cluster2","nonPJI_cluster1","nonPJI_cluster2"))

UMAPdata_timecourse <- na.omit(UMAPdata_all[,c("ID","UMAP1","UMAP2","PJI","PJI_cluster","Time_original")])
UMAPdata_timecourse$PJI <- factor(UMAPdata_timecourse$PJI,levels=c("1","0"))

p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_timecourse) +
  geom_point(aes(fill=PJI),shape=21,alpha=0.8,size=3) +
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  #scale_colour_brewer(palette = "Set1") + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  facet_grid(PJI ~ Time_original)+
  scale_fill_manual(values=c("red","gray"))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)
ggsave(file = "YCU_PJI_postope_SRFdist_UMAP_diagnosis.pdf", plot = p,width=12,height=12)


p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_timecourse) +
  geom_path(aes(colour=PJI,group=ID),linetype=2,arrow = arrow(type = "closed"))+
  geom_point(aes(fill=PJI),shape=21,alpha=0.8,size=3) +
  geom_polygon(data=ellipse,aes(group=colour),color="black",fill=NA,linetype=2)+
  #scale_colour_brewer(palette = "Set1") + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  # facet_grid(PJI_cluster ~ Time_original)+
  scale_fill_manual(values=c("gray","red"))+
  scale_colour_manual(values=c("gray","red"))+
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_bw(base_size = 16)+
  theme(aspect.ratio=1)
ggsave(file = "YCU_PJI_postope_SRFdist_UMAP_diagnosis_path.pdf", plot = p,width=12,height=12)


continuous_vars <- colnames(original_data)[6:ncol(original_data)]

for (item in continuous_vars){
  UMAPdata_ex <- na.omit(data.frame(UMAPdata_postope[,c("UMAP1","UMAP2")],UMAPdata_postope[,item,drop=F]))
  colnames(UMAPdata_ex) <- c("UMAP1","UMAP2","value")

  u_lim <- quantile(UMAPdata_ex$value,0.95)
  UMAPdata_ex$value <- sapply(UMAPdata_ex$value,u_lim=u_lim,cutoff)

  p <- ggplot(aes(x=UMAP1, y=UMAP2),data=UMAPdata_ex) +
  geom_point(aes(color=value),shape=16,alpha=0.8,size=3) +
  scale_colour_distiller(palette = "Spectral") +
  # scale_colour_manual(values = c("gray","red")) + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  #scale_fill_manual(values=rev(cols[c(1,2)]))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)
  ggsave(file = paste0("YCU_PJI_postope_SRFdist_UMAP_",item,".pdf"), plot = p,width=7,height=6)
}

UMAPdata_merged <- data.frame(rbind(UMAPdata_preope,UMAPdata_postope))
UMAPdata_merged$Time_original <- factor(UMAPdata_merged$Time_original, levels=c("1m_ago","3m_later","6m_later","12m_later"))


for (item in continuous_vars){
  UMAPdata_ex <- na.omit(data.frame(UMAPdata_merged[,c("UMAP1","UMAP2")],UMAPdata_merged[,item,drop=F]))
  colnames(UMAPdata_ex) <- c("UMAP1","UMAP2","value")

  u_lim <- quantile(UMAPdata_ex$value,0.95)
  UMAPdata_ex$value <- sapply(UMAPdata_ex$value,u_lim=u_lim,cutoff)

  p <- ggplot(aes(x=UMAP1, y=UMAP2),data=UMAPdata_ex) +
  geom_point(aes(color=value),shape=16,alpha=0.8,size=3) +
  scale_colour_distiller(palette = "Spectral") +
  # scale_colour_manual(values = c("gray","red")) + 
  #scale_shape_manual(values=c(1,17))+
  #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
  #              geom="polygon", level=0.95, alpha=0.2) +
  #scale_fill_manual(values=rev(cols[c(1,2)]))+
  # scale_colour_manual(values=cols) +
  scale_x_continuous(limits = c(min(UMAPdata_preope$UMAP1)*1.1 - max(UMAPdata_preope$UMAP1)*0.15,max(UMAPdata_preope$UMAP1)*1.15-min(UMAPdata_preope$UMAP1)*0.1)) +
  scale_y_continuous(limits = c(min(UMAPdata_preope$UMAP2)*1.1 - max(UMAPdata_preope$UMAP2)*0.15,max(UMAPdata_preope$UMAP2)*1.15-min(UMAPdata_preope$UMAP2)*0.1)) +
  theme_classic(base_size = 16)
  ggsave(file = paste0("YCU_PJI_merged_SRFdist_UMAP_",item,".pdf"), plot = p,width=7,height=6)
}

rownames(UMAPdata_merged) <- c(rownames(original_data_train),rownames(original_data_test))

antibacteria_bacteria_PJI <- read.table("pji_koukinyaku_plus-bg1_bacteria_PJI.txt",sep="\t",header=T,row.names=1)
antibacteria_bacteria_list <- as.vector(colnames(antibacteria_bacteria_PJI)[c(9:23)])

postope_timeseries <- c("3m_later","6m_later","12m_later")

for (i in 1:3){

  UMAPdata_time <- UMAPdata_merged[UMAPdata_merged$Time_original %in% c("1m_ago",postope_timeseries[i]),]
  col <- cols2[c(2:4)][i]

  p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_time[UMAPdata_time$PJI==1,]) +
    geom_path(aes(group=ID),color=col,linetype=2,alpha=0.8,,arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches")))+
    geom_point(aes(color=Time_original),shape=15,alpha=0.8,size=3,stroke=1) +
    scale_shape_manual(values=c(1,15))+
    scale_colour_manual(values=c("gray",col)) + 
    #scale_shape_manual(values=c(1,17))+
    #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
    #              geom="polygon", level=0.95, alpha=0.2) +
    #scale_colour_manual(values=cols[c(2,3)])+
    # scale_colour_manual(values=cols) +
    scale_x_continuous(limits = c(min(UMAPdata_merged$UMAP1)*1.1 - max(UMAPdata_merged$UMAP1)*0.1,max(UMAPdata_merged$UMAP1)*1.1-min(UMAPdata_merged$UMAP1)*0.1)) +
    scale_y_continuous(limits = c(min(UMAPdata_merged$UMAP2)*1.1 - max(UMAPdata_merged$UMAP2)*0.1,max(UMAPdata_merged$UMAP2)*1.1-min(UMAPdata_merged$UMAP2)*0.1)) +
    theme_classic(base_size = 16)
  ggsave(file = paste0("YCU_PJI_SRFdist_UMAP_PJI_postope_",postope_timeseries[i],".pdf"), plot = p,width=7,height=6)

  UMAPdata_time_PJI <- UMAPdata_time[UMAPdata_time$PJI==1,]

  for (ab in antibacteria_bacteria_list){

    pat_list <- as.vector(rownames(na.omit(antibacteria_bacteria_PJI[,ab,drop=F])))

    UMAPdata_time_PJI_ab <- UMAPdata_time_PJI[UMAPdata_time_PJI$ID %in% pat_list,]

    p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_time_PJI_ab) +
      geom_path(aes(group=ID),color=col,linetype=2,alpha=0.8,,arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches")))+
      geom_point(aes(color=Time_original),shape=15,alpha=0.8,size=3,stroke=1) +
      scale_shape_manual(values=c(1,15))+
      scale_colour_manual(values=c("gray",col)) + 
      #scale_shape_manual(values=c(1,17))+
      #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
      #              geom="polygon", level=0.95, alpha=0.2) +
      #scale_colour_manual(values=cols[c(2,3)])+
      # scale_colour_manual(values=cols) +
      scale_x_continuous(limits = c(min(UMAPdata_merged$UMAP1)*1.1 - max(UMAPdata_merged$UMAP1)*0.1,max(UMAPdata_merged$UMAP1)*1.1-min(UMAPdata_merged$UMAP1)*0.1)) +
      scale_y_continuous(limits = c(min(UMAPdata_merged$UMAP2)*1.1 - max(UMAPdata_merged$UMAP2)*0.1,max(UMAPdata_merged$UMAP2)*1.1-min(UMAPdata_merged$UMAP2)*0.1)) +
      theme_classic(base_size = 16)
    ggsave(file = paste0("YCU_PJI_SRFdist_UMAP_PJI_postope_",postope_timeseries[i],"_",ab,".pdf"), plot = p,width=7,height=6)
  }

  p <- ggplot(aes(x=UMAP1,y=UMAP2),data=UMAPdata_time[UMAPdata_time$PJI==0,]) +
    geom_path(aes(group=ID),color=col,linetype=2,arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches")))+
    geom_point(aes(color=Time_original),shape=1,alpha=0.8,size=3,stroke=1) +
    scale_shape_manual(values=c(1,15))+
    scale_colour_manual(values=c("gray",col)) + 
    #scale_shape_manual(values=c(1,17))+
    #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
    #              geom="polygon", level=0.95, alpha=0.2) +
    #scale_colour_manual(values=cols[c(2,3)])+
    # scale_colour_manual(values=cols) +
    scale_x_continuous(limits = c(min(UMAPdata_merged$UMAP1)*1.1 - max(UMAPdata_merged$UMAP1)*0.1,max(UMAPdata_merged$UMAP1)*1.1-min(UMAPdata_merged$UMAP1)*0.1)) +
    scale_y_continuous(limits = c(min(UMAPdata_merged$UMAP2)*1.1 - max(UMAPdata_merged$UMAP2)*0.1,max(UMAPdata_merged$UMAP2)*1.1-min(UMAPdata_merged$UMAP2)*0.1)) +
    theme_classic(base_size = 16)
  ggsave(file = paste0("YCU_PJI_SRFdist_UMAP_nonPJI_postope_",postope_timeseries[i],".pdf"), plot = p,width=7,height=6)
}


for (pID in unique(UMAPdata_merged$ID)){
  UMAPdata_ind <- UMAPdata_merged[UMAPdata_merged$ID==pID,]
  UMAPdata_ind$Time_original <- factor(UMAPdata_ind$Time_original, levels=c("1m_ago","3m_later","6m_later","12m_later"))
  UMAPdata_ind <- UMAPdata_ind[order(UMAPdata_ind$Time_original),]

  if (UMAPdata_ind$PJI[[1]] == "0"){
    # col = cols[1]
    folder = "nonPJI_patients"
    shape = 1
  } else {
    # col = cols[2]
    folder = "PJI_patients"
    shape = 15
  }

  p <- ggplot() +
    geom_point(aes(x=UMAP1,y=UMAP2,shape=PJI),data=UMAPdata_preope,color="gray",size=3)+
    geom_point(aes(x=UMAP1,y=UMAP2,color=Time_original),data=UMAPdata_ind,shape=shape,size=3,stroke=1) +
    geom_path(aes(x=UMAP1,y=UMAP2),data=UMAPdata_ind,color="gray40", linetype=2,arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) +
    scale_shape_manual(values=c(1,15))+
    # geom_text_repel(aes(label = date)) +
    # scale_shape_manual(values=c(16,1,2,7),drop = FALSE)+
    #scale_colour_brewer(palette = "Set1") + 
    #scale_shape_manual(values=c(1,17))+
    #stat_ellipse(aes(x=MDS1,y=MDS2,fill=cluster),
    #              geom="polygon", level=0.95, alpha=0.2) +
    scale_colour_manual(values=cols2,drop=F)+
    # scale_colour_manual(values=cols) +
    scale_x_continuous(limits = c(min(UMAPdata_merged$UMAP1)*1.1 - max(UMAPdata_merged$UMAP1)*0.1,max(UMAPdata_merged$UMAP1)*1.1-min(UMAPdata_merged$UMAP1)*0.1)) +
    scale_y_continuous(limits = c(min(UMAPdata_merged$UMAP2)*1.1 - max(UMAPdata_merged$UMAP2)*0.1,max(UMAPdata_merged$UMAP2)*1.1-min(UMAPdata_merged$UMAP2)*0.1)) +
    theme_classic(base_size = 16)
  ggsave(file = paste0("./",folder,"/YCU_PJI_",pID,"_SRFdist_UMAP.pdf"), plot = p,width=7,height=6)
}
