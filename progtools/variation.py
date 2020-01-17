'''
VARIATION.PY
	This file contains methods for assessing variation among NPs,
	including generation of a normalized similarity matrix and a PCA
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import decomposition
import pandas as pd
import numpy as np
import seaborn as sns

# SIMILARITY: takes as input a segmentation table
#	generates a normalized similarity matrix between NPs
#	normalization process:
#		1. calculate a similarity (summed AND operator) matrix
#		2. compute total number of 1s for each row
#		3. given a pairwise comparison, the normalized value is the similarity divided by the lowest total num 1s
def similarity(dfseg,outdir):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	st = st.T # make NPs the index (to compute similarity matrix of windows, use original dfseg)
	simmat = (st.values & st.values[:, None]).sum(2) # summed AND operation to make unnormalized similarity
	rowsums = np.sum(st.values,axis=1) # sum up windows for each row
	normsim = []
	for x in range(0,len(simmat)): # iterate over pairwise NPs and calculate the normalized value
		newr = []
		for y in range(0,len(simmat)):
			minval = min((rowsums[x],rowsums[y])) # get lowest total num 1s
			if(minval == 0):
				newr.append(0) # edge case: avoid div by 0
			else:
				newr.append(simmat[x][y]/minval) # similarity divided by the lowest total num 1s gives normalized
		normsim.append(newr)
	simdf = pd.DataFrame(normsim,columns=st.index,index=st.index) # np array
	simdf.to_csv(outdir + "NP-similarity-matrix.csv")
	print("-- Finished generating normalized similarity matrix, saved to",outdir + "NP-similarity-matrix.csv")
	fig, ax = plt.subplots(figsize=(10,10),facecolor='white')
	sns.clustermap(data=simdf,metric="jaccard",method="complete")
	plt.title("Nuclear Profile Similarity Matrix")
	plt.savefig(outdir + "NP-similarity-heatmap.png")
	plt.clf()
	print("-- Finished generating heatmap of normalized similarity matrix, saved to",outdir + "NP-similarity-heatmap.png")

# SIMILARITY: takes as input a segmentation table
#	generates a normalized similarity matrix between NPs
def pca(dfseg,outdir):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)].T # drop NPS which do not cover anything in GRI and set rows to NPs
	arr = st.values # get simple np array of values
	numcols = len(arr[0]) - 1
	pca = decomposition.PCA(n_components=2).fit_transform(arr) # create model and fit to the data
	fig, ax = plt.subplots(figsize=(10,10))
	w = [] # first principal component
	r = [] # second principal component
	for row in pca:
		w.append(row[0])
		r.append(row[1])  
	ax.scatter(w, r) # plot pca to 2d scatterplot
	ax.set_xlabel('First Principal Component')
	ax.set_ylabel('Second Principal Component')
	plt.title("2d PCA of Nuclear Profiles in GRI")
	plt.savefig(outdir + "NP-PCA.png")
	plt.clf()
	print("-- Finished generating PCA of nuclear profiles, saved to",outdir + "NP-PCA.png")