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
UMAPdata_all$cluster <- clu

UMAPdata_preope$PJI <- as.factor(UMAPdata_preope$PJI)

UMAPdata_all$PJI <- as.factor(UMAPdata_all$PJI)
UMAPdata_all$Time_original <- factor(UMAPdata_all$Time_original, levels=c("1m_ago","3m_later","6m_later","12m_later"))

UMAPdata_all_wider <- UMAPdata_all[,c("ID","PJI","Time_original","cluster")] %>% pivot_wider(id_cols=c(ID,PJI),names_from=Time_original,values_from=cluster)
UMAPdata_all_wider$`1m_ago` <- UMAPdata_all_wider$`1m_ago` + 2
UMAPdata_all_wider$`3m_later` <- UMAPdata_all_wider$`3m_later` + 4
UMAPdata_all_wider$`6m_later` <- UMAPdata_all_wider$`6m_later` + 6
UMAPdata_all_wider$`12m_later` <- UMAPdata_all_wider$`12m_later` + 8


edge_1mago <- UMAPdata_all_wider %>% dplyr::count(PJI,`1m_ago`) %>% na.omit() %>% as_tbl_graph()

edge_1mago_to_3mlater <- UMAPdata_all_wider %>% dplyr::count(`1m_ago`, `3m_later`) %>% na.omit() %>% as_tbl_graph()

edge_3mlater_to_6mlater <- UMAPdata_all_wider %>% dplyr::count(`3m_later`, `6m_later`) %>% na.omit() %>% as_tbl_graph()

edge_6mlater_to_12mlater <- UMAPdata_all_wider %>% dplyr::count(`6m_later`, `12m_later`) %>% na.omit() %>% as_tbl_graph()

g <- edge_1mago %>% graph_join(edge_1mago_to_3mlater) %>% graph_join(edge_3mlater_to_6mlater) %>% graph_join(edge_6mlater_to_12mlater)

edge <- g %>% activate(edges) %>% as_tibble() %>% mutate(from=from-1,to=to-1)
node <- g %>% activate(nodes) %>% as_tibble()

p2 <- sankeyNetwork(
  Links = edge,
  Nodes = node,
  Source = "from",
  Target = "to",
  Value = "n",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30,
  iterations = 0
)

saveNetwork(p2, file="PJI_state_transition_sankeyplot.html")



UMAPdata_PJI <- UMAPdata_all[UMAPdata_all$PJI=="1",]
UMAPdata_PJI_wider <- UMAPdata_PJI[,c("ID","PJI","Time_original","cluster")] %>% pivot_wider(id_cols=c(ID,PJI),names_from=Time_original,values_from=cluster)
UMAPdata_PJI_wider$`1m_ago` <- UMAPdata_PJI_wider$`1m_ago` + 2
UMAPdata_PJI_wider$`3m_later` <- UMAPdata_PJI_wider$`3m_later` + 4
UMAPdata_PJI_wider$`6m_later` <- UMAPdata_PJI_wider$`6m_later` + 6
UMAPdata_PJI_wider$`12m_later` <- UMAPdata_PJI_wider$`12m_later` + 8

PJI_1mago_edge <- UMAPdata_PJI_wider %>% 
  dplyr::count(PJI,`1m_ago`) %>% na.omit() %>% 
  as_tbl_graph()

PJI_1mago_to_3mlater_edge <- UMAPdata_PJI_wider %>% 
  dplyr::count(`1m_ago`, `3m_later`) %>% na.omit() %>% 
  as_tbl_graph()

PJI_3mlater_to_6mlater_edge <- UMAPdata_PJI_wider %>% 
  dplyr::count(`3m_later`, `6m_later`) %>% na.omit() %>% 
  as_tbl_graph()

PJI_6mlater_to_12mlater_edge <- UMAPdata_PJI_wider %>% 
  dplyr::count(`6m_later`, `12m_later`) %>% na.omit() %>% 
  as_tbl_graph()

g <- PJI_1mago_edge %>% graph_join(PJI_1mago_to_3mlater_edge) %>% graph_join(PJI_3mlater_to_6mlater_edge) %>% graph_join(PJI_6mlater_to_12mlater_edge)

edge <- g %>% activate(edges) %>% as_tibble() %>% mutate(from=from-1,to=to-1)
node <- g %>% activate(nodes) %>% as_tibble()

p2 <- sankeyNetwork(
  Links = edge,
  Nodes = node,
  Source = "from",
  Target = "to",
  Value = "n",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30,
  iterations = 0
)

saveNetwork(p2, file="PJI_state_transition_sankeyplot_PJI.html")


UMAPdata_nonPJI <- UMAPdata_all[UMAPdata_all$PJI=="0",]
UMAPdata_nonPJI_wider <- UMAPdata_nonPJI[,c("ID","PJI","Time_original","cluster")] %>% pivot_wider(id_cols=c(ID,PJI),names_from=Time_original,values_from=cluster)
UMAPdata_nonPJI_wider$`1m_ago` <- UMAPdata_nonPJI_wider$`1m_ago` + 2
UMAPdata_nonPJI_wider$`3m_later` <- UMAPdata_nonPJI_wider$`3m_later` + 4
UMAPdata_nonPJI_wider$`6m_later` <- UMAPdata_nonPJI_wider$`6m_later` + 6
UMAPdata_nonPJI_wider$`12m_later` <- UMAPdata_nonPJI_wider$`12m_later` + 8

nonPJI_1mago_edge <- UMAPdata_nonPJI_wider %>% 
  dplyr::count(PJI,`1m_ago`) %>% na.omit() %>% 
  as_tbl_graph()

nonPJI_1mago_to_3mlater_edge <- UMAPdata_nonPJI_wider %>% 
  dplyr::count(`1m_ago`, `3m_later`) %>% na.omit() %>% 
  as_tbl_graph()

nonPJI_3mlater_to_6mlater_edge <- UMAPdata_nonPJI_wider %>% 
  dplyr::count(`3m_later`, `6m_later`) %>% na.omit() %>% 
  as_tbl_graph()

nonPJI_6mlater_to_12mlater_edge <- UMAPdata_nonPJI_wider %>% 
  dplyr::count(`6m_later`, `12m_later`) %>% na.omit() %>% 
  as_tbl_graph()

g <- nonPJI_1mago_edge %>% graph_join(nonPJI_1mago_to_3mlater_edge) %>% graph_join(nonPJI_3mlater_to_6mlater_edge) %>% graph_join(nonPJI_6mlater_to_12mlater_edge)

edge <- g %>% activate(edges) %>% as_tibble() %>% mutate(from=from-1,to=to-1)
node <- g %>% activate(nodes) %>% as_tibble()

p2 <- sankeyNetwork(
  Links = edge,
  Nodes = node,
  Source = "from",
  Target = "to",
  Value = "n",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30,
  iterations = 0
)

saveNetwork(p2, file="PJI_state_transition_sankeyplot_nonPJI.html")


