library(markovchain)


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

PJIpatients <- unique(UMAPdata_all[UMAPdata_all$PJI==1,"ID"])
nonPJIpatients <- unique(UMAPdata_all[UMAPdata_all$PJI==0,"ID"])

PJIsequences <- list()
nonPJIsequences <- list()

time_point <- c("1m_ago","3m_later","6m_later","9m_later","12m_later")

for (ID in PJIpatients){
	ind_data <- UMAPdata_all[UMAPdata_all$ID==ID,]
	tmp_seq <- c()
	for (t in time_point){
		if (t %in% ind_data$Time_original){
			tmp_seq <- c(tmp_seq,ind_data[ind_data$Time_original==t,"cluster"])
		} else {
			tmp_seq <- c(tmp_seq,NA)
		}
	}
	PJIsequences[ID] <- list(tmp_seq)
}

for (ID in nonPJIpatients){
	ind_data <- UMAPdata_all[UMAPdata_all$ID==ID,]
	tmp_seq <- c()
	for (t in time_point){
		if (t %in% ind_data$Time_original){
			tmp_seq <- c(tmp_seq,ind_data[ind_data$Time_original==t,"cluster"])
		} else {
			tmp_seq <- c(tmp_seq,NA)
		}
	}
	nonPJIsequences[ID] <- list(tmp_seq)
}

PJI_mcFitMle <- markovchainFit(PJIsequences,method="mle")
nonPJI_mcFitMle <- markovchainFit(nonPJIsequences,method="mle")


