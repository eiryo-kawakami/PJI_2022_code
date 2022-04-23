import pandas as pd
import numpy as np
import SRFdist as srfd
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
	res = srfd.SRFdist(datRF_sep_train, no_tree=100, no_rep=10, imp=True, max_depth=5, njobs=6 )
	RFimp = pd.concat([RFimp, res["RFimportance"]], axis=1)
	for j in range(len(res["models"])):
		with open("./SRFdist_models/YCU_PJI_preope_imp"+str(i+1)+"_SRFdist_rep"+str(j+1)+".sav", 'wb') as f:
			cPickle.dump(res["models"][j], f)
		test_prob = pd.DataFrame(res["models"][j].predict_proba(datRF_sep_test)).iloc[:,1]
		print(test_prob)
		dat_info_test["prob_"+str(i)+"_"+str(j)] = list(test_prob)
		leaves = res["models"][j].apply(datRF_sep.loc[:,list(exp_vars)])
		leaves_df = pd.concat([leaves_df, pd.DataFrame(leaves)], axis=1)

# leaves_df.index = datRF.index
print(leaves_df)

leaves_df_train = leaves_df.iloc[train_list,:]
leaves_df_test = leaves_df.iloc[test_list,:]

leaves_df_train.to_csv("YCU_PJI_preope_train_SRFdist_python_leaves.txt",sep="\t")
leaves_df_test.to_csv("YCU_PJI_preope_test_SRFdist_python_leaves.txt",sep="\t")

dat_info_test.to_csv("YCU_PJI_preope_SRFprob.txt",sep="\t")
RFimp.to_csv("YCU_PJI_preope_SRFimportance.txt",sep="\t")
