import random
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.utils import check_random_state
import joblib

def SRFdist(datRF, no_tree=4000, no_rep=10, imp=True, max_depth=None ,proxConver=False, njobs=1 ):

	def _generate_sample_indices(random_state, n_samples):
		"""Private function used to _parallel_build_trees function."""
		random_instance = check_random_state(random_state)
		sample_indices = random_instance.randint(0, n_samples, n_samples)

		return sample_indices

	def _generate_unsampled_indices(random_state, n_samples):
		"""Private function used to forest._set_oob_score function."""
		sample_indices = _generate_sample_indices(random_state, n_samples)
		sample_counts = np.bincount(sample_indices, minlength=n_samples)
		unsampled_mask = sample_counts == 0
		indices_range = np.arange(n_samples)
		unsampled_indices = indices_range[unsampled_mask]

		return unsampled_indices

	def calc_OOB(RFmodel,tree_id,nrow1):
		OOB_samples = _generate_sample_indices(RF1.estimators_[0].random_state,nrow1)
		is_OOB = np.array([ 1 if k in OOB_samples else 0 for k in np.array(range(nrow1))])

		return is_OOB

	def cleandist(x):
		x = np.where(x <= 0.0, 0.0000000001, x)
		return x

	res = {}

	nrow1 = len(datRF)
	ncol1 = len(datRF.columns)
	RFprox = np.zeros((nrow1, nrow1))
	RFimportance = pd.DataFrame([])

	res["RFimportance"] = pd.DataFrame([])
	res["leaves_df"] = pd.DataFrame([])

	y = list(datRF.iloc[:,0])
	# print(datRF)
	datRFx = datRF.iloc[:,1:]
	rep1 = [666] * nrow1
	res["prob"] = np.array([0.0] * nrow1)
	res["models"] = []

	for i in range(no_rep):
		print(i)
		index1 = np.array(range(nrow1))
		random.shuffle(index1)
		yy = []
		for j in range(nrow1):
			rep1[index1[j]] = j
			yy.append(y[index1[j]])

		datRFsyn = datRFx.iloc[index1,]

		clf = RandomForestClassifier(n_estimators=no_tree,criterion="gini",max_features="auto",max_depth=max_depth, random_state=0, n_jobs=njobs)
		clf_isotonic = CalibratedClassifierCV(clf, method='isotonic')
		RF1 = clf_isotonic.fit(datRFsyn, yy)
		res["models"].append(clf)
		leaves = RF1.apply(datRFsyn)[rep1,:]
		res["leaves_df"] = pd.concat([res["leaves_df"], pd.DataFrame(leaves)], axis=1)
		res["prob"] += np.array(pd.DataFrame(clf.predict_proba(datRFsyn)).iloc[:,1][rep1])

		if imp:
			res["RFimportance"] = pd.concat([res["RFimportance"], pd.DataFrame(RF1.feature_importances_)], axis=1)

	# distRF = cleandist(np.sqrt(np.ones((nrow1,nrow1))-RFprox/no_rep))
	# distRF = pd.DataFrame(distRF)
	# distRF.index = datRF.index
	# distRF.columns = datRF.index

	res["prob"] = res["prob"] / no_rep

	if imp:
		res["RFimportance"].columns = [ "rep_"+str(i+1) for i in range(no_rep)]
		res["RFimportance"].index = datRFx.columns

	return res
