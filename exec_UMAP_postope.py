import umap
import numpy as np
import pandas as pd
import _pickle as cPickle


df = pd.read_csv("YCU_PJI_postope_SRFdist_python_leaves.txt",sep="\t",header=0,index_col=0)

with open("YCU_PJI_preope_SRFdist_UMAP.sav", 'rb') as f:
	trans = cPickle.load(f)
	embedding = pd.DataFrame(trans.transform(df))
embedding.to_csv("YCU_PJI_postope_SRFdist_UMAP_embedding.txt",sep="\t")
