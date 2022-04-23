import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import _pickle as cPickle

n_imp = 10
n_rep = 10

dat = pd.read_csv("PJIdata_MICEimp_postope.txt",sep="\t",header=0)
exp_vars = dat.columns[5:len(dat.columns)]
datRF = dat.loc[:,exp_vars]

sample_num = int(len(datRF)/10)

leaves_df = pd.DataFrame([])

for i in range(n_imp):
	datRF_sep = datRF.iloc[(i*sample_num):((i+1)*sample_num),:]
	for j in range(n_rep):
		with open("./SRFdist_models/YCU_PJI_preope_imp"+str(i+1)+"_SRFdist_rep"+str(j+1)+".sav", 'rb') as f:
			model = cPickle.load(f)
			leaves = model.apply(datRF_sep)
			leaves_df = pd.concat([leaves_df, pd.DataFrame(leaves)], axis=1)

leaves_df.to_csv("YCU_PJI_postope_SRFdist_python_leaves.txt",sep="\t")
# RFimportance1.to_csv("Chiba_PK_all_RFimportance_python.txt",sep="\t")
