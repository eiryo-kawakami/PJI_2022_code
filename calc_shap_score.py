import shap
from sklearn.ensemble import RandomForestClassifier
import _pickle as cPickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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

datRF_sep = datRF.iloc[0:sample_num,:]
datRF_sep_train = datRF_sep.iloc[train_list,:]
X_train = datRF_sep_train.loc[:,list(exp_vars)]

shap.initjs()

# def shap_summaryplot(model,X):
# 	explainer = shap.TreeExplainer(model)
# 	shap_values = explainer.shap_values(X)

# 	return shap.summary_plot(shap_values[1], X, show=False)

with open("./SRFdist_models/YCU_PJI_preope_imp1_SRFdist_rep1.sav", 'rb') as f:
	model = cPickle.load(f)

	explainer = shap.TreeExplainer(model)
	shap_values = explainer.shap_values(X_train)

	plt.clf()
	shap.summary_plot(shap_values[1], X_train, show=False)
	plt.savefig("YCU_PJI_preope_imp1_SRFdist_rep1_shap_summaryplot.pdf")

	plt.clf()
	shap.dependence_plot('CRP', shap_values[1], X_train, interaction_index='CRP', show=False)
	plt.savefig("YCU_PJI_preope_imp1_SRFdist_rep1_shap_dependence_plot_CRP.pdf")

	shap_interaction_values = explainer.shap_interaction_values(X_train)

	plt.clf()
	shap.dependence_plot(('CRP','Total.protein'), shap_interaction_values[1], X_train, show=False)
	plt.savefig("YCU_PJI_preope_imp1_SRFdist_rep1_shap_dependence_plot_CRPvsTP.pdf")


	
