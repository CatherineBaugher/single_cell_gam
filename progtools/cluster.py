'''
PREP.PY
	This file contains clustering steps,
	including generating a clustered heatmap and correlation matrix,
	and finally calculating radial position and compaction
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn import decomposition
from sklearn.cluster import KMeans
import seaborn as sns
import pandas as pd
from random import random
import numpy as np

# COMPACTION: takes as input a list of list containing NP clusters
#	saves a boxplot of compaction (how many NPs in a cluster each window of the GRI contains)
#	writes avgs of each cluster to standard output
def compaction(dfseg,npclusts,outf):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	collection = [] # store compaction calculations
	for x in range(0,len(npclusts)):
		selected = st[npclusts[x]] # constrain seg table to only nps in the cluster
		totnps = len(npclusts[x])
		out1 = selected.sum(axis=1) # calculate sum
		out1.to_csv(outf + "cluster" + str(x + 1) + "-compaction-unnormalized.csv")
		out2 = out1 / totnps # normalize sum by total num NPs in the cluster
		#out2 = out2[out2 >= .25] # OPTIONAL? -- get only windows with >= .25
		out2.to_csv(outf + "cluster" + str(x + 1) + "-compaction-normalized.csv")
		collection.append(out2.values)
		print("-- Cluster",str(x+1),"has compaction average",sum(out2.values) / len(out2))
	print("-- Finished calculating cluster compactions, unnormalized and normalized values saved to files in form",outf + "clusterX-compaction-normalized.csv","and",outf + "clusterX-compaction-unnormalized.csv")
	fig, ax = plt.subplots(figsize=(10,10),facecolor='white')
	ax.boxplot(collection,labels=["Cluster " + str(x+1) for x in range(0,len(npclusts))])
	plt.title("Boxplot of Clustered Compactions Across GRI")
	plt.savefig(outf + "cluster-compaction-boxplot.png")
	print("-- Finished generating boxplot of cluster compactions, saved to",outf + "cluster-compaction-boxplot.png")

# KMEANSCLUST: takes as input a segmentation table and optional threshold for number of clusters
#	saves a PCA of NPs clustering into clustparam number of groups
def kmeansclust(dfseg,outf,clustparam=3):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)].T # drop NPS which do not cover anything in GRI and set rows to NPs
	# generate PCA
	arr = st.values # get simple np array of values
	numcols = len(arr[0]) - 1
	pca = decomposition.PCA(n_components=2).fit_transform(arr) # create model and fit to the data
	fig, ax = plt.subplots(figsize=(10,10))
	# calculate k-means clusters
	pts = np.array([[row[0],row[1]] for row in pca])
	kmeansclust = KMeans(n_clusters=clustparam, random_state=0).fit(pts)
	kmlabels = kmeansclust.labels_
	centroids = kmeansclust.cluster_centers_
	clustlabs = pd.DataFrame(columns=["cluster"],index=st.index)
	clustlabs["cluster"] = kmlabels
	for i in range(clustparam):
		thisclust = pts[np.where(kmlabels==i)] # find nps with current cluster label
		plt.plot(thisclust[:,0],thisclust[:,1],'o',label="Cluster "+str(i+1)) # plot the data observations
		plt.plot(centroids[i,0],centroids[i,1],'kx') # plot centroid
    # create figure
	ax.set_xlabel('First Principal Component')
	ax.set_ylabel('Second Principal Component')
	plt.title("K-means Clustering of Nuclear Profiles - 2d PCA")
	plt.savefig(outf + "cluster-kmeans-pca.png")
	clustlabs.to_csv(outf + "cluster-kmeans-raw.csv")
	print("-- Finished generating K-means clustering, saved to",outf + "cluster-kmeans-raw.csv")
	print("-- Finished generating PCA of K-means clustering, saved to",outf + "cluster-kmeans-pca.png")
	listofclusts = [[] for x in range(0,clustparam)]
	for ind,row in clustlabs.iterrows():
		listofclusts[row["cluster"]].append(ind)
	return listofclusts

# HEATMAPCLUST: takes as input a segmentation table and optional threshold for max number of clusters
#	saves a clustermap of NPs using jaccard distance
#	precomputed clusters may also be fed to the program to color labels
def heatmapclust(dfseg,outf,clustlabs,ctype="single"):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	st = st.T # make NPs be the rows
	stdist = distance.pdist(st,metric="jaccard")
	dendro = hierarchy.linkage(stdist,method=ctype)
	numclust = len(clustlabs) # take the number of clusters from the list of lists
	mycolors = [] # get set of N=numclust number of colors for labelling clustermap
	r = int(random() * 256) # initialize modulus values for red, green, blue
	g = int(random() * 256)
	b = int(random() * 256)
	for i in range(1,numclust+1):
		r = (r + (256/numclust)) % 256
		g = (g + (256/numclust)) % 256
		b = (b + (256/numclust)) % 256
		mycolors.append((r/256,g/256,b/256)) # constrain values 0-1
	classinfo = {} # make dictionary of NPs and their respective clusters
	for c in range(0,len(clustlabs)):
		for np in clustlabs[c]:
			classinfo[np] = mycolors[c]
	figname = "cluster-heatmap-"+ctype+".png"
	fig, ax = plt.subplots(figsize=(10,10),facecolor="white")
	myclust = sns.clustermap(st,row_linkage=dendro,col_cluster=False) # use linkage to build clustermap
	for lab in myclust.ax_heatmap.axes.get_yticklabels():
		npname = lab.get_text()
		lab.set_color(classinfo[npname])
	plt.gcf().subplots_adjust(bottom=0.20)
	plt.suptitle("Clustered Heatmap of NPs")
	plt.savefig(outf + figname)
	print("-- Finished generating clustered heatmap, saved to",outf + figname)