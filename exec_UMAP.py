import umap
import numpy as np
import pandas as pd
import _pickle as cPickle


df_train = pd.read_csv("YCU_PJI_preope_train_SRFdist_python_leaves.txt",sep="\t",header=0,index_col=0)
train_num = len(df_train.index)

df_test = pd.read_csv("YCU_PJI_preope_test_SRFdist_python_leaves.txt",sep="\t",header=0,index_col=0)
test_num = len(df_test.index)

df_postope = pd.read_csv("YCU_PJI_postope_SRFdist_python_leaves.txt",sep="\t",header=0,index_col=0)
postope_num = len(df_postope.index)

df = pd.concat([df_train, df_test, df_postope], axis=0)

reducer = umap.UMAP(metric="hamming",min_dist=0.001)
trans = reducer.fit(X=df.values)
with open("YCU_PJI_SRFdist_UMAP.sav", 'wb') as f:
	cPickle.dump(trans, f)
embedding = pd.DataFrame(trans.embedding_)

embedding_train = embedding.iloc[0:train_num,:]
embedding_test = embedding.iloc[train_num:(train_num+test_num),:]
embedding_postope = embedding.iloc[(train_num+test_num):(train_num+test_num+postope_num),:]

embedding_train.to_csv("YCU_PJI_preope_train_SRFdist_UMAP_embedding.txt",sep="\t")
embedding_test.to_csv("YCU_PJI_preope_test_SRFdist_UMAP_embedding.txt",sep="\t")
embedding_postope.to_csv("YCU_PJI_postope_SRFdist_UMAP_embedding.txt",sep="\t")


