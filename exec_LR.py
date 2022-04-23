import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
import _pickle as cPickle

n_imp = 10

dat = pd.read_csv("PJIdata_MICEimp_preope.txt",sep="\t",header=0)
exp_vars = dat.columns[5:len(dat.columns)]
datRF = dat.loc[:,["PJI"]+list(exp_vars)]

sample_num = int(len(datRF)/10)
train_list = [a for a in range(sample_num) if a % 3 != 2]
test_list = [a for a in range(sample_num) if a % 3 == 2]

leaves_df = pd.DataFrame([])
dat_info_test = dat.iloc[test_list,[0,1,2]]

RFimp = pd.DataFrame([])

for i in range(n_imp):
	datRF_sep = datRF.iloc[(i*sample_num):((i+1)*sample_num),:]
	datRF_sep_train = datRF_sep.iloc[train_list,:]
	datRF_sep_test = datRF_sep.iloc[test_list,:]
	datRF_sep_test = datRF_sep_test.loc[:,list(exp_vars)]
	res = LogisticRegression(penalty="none", n_jobs=6, max_iter=1000).fit(datRF_sep_train.iloc[:,1:],datRF_sep_train.iloc[:,0])
	# print(res.predict_proba(datRF_sep_test))
	test_prob = pd.DataFrame(res.predict_proba(datRF_sep_test)).iloc[:,1]
	dat_info_test["prob_"+str(i)] = list(test_prob)

dat_info_test.to_csv("YCU_PJI_preope_LRprob.txt",sep="\t")