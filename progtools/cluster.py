'''
PREP.PY
	This file contains clustering steps,
	including generating a clustered heatmap and correlation matrix,
	and finally calculating radial position and compaction
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# COMPACTION: takes as input a list of list containing NP clusters
#	saves a boxplot of compaction (how many NPs in a cluster each window of the GRI contains)
#	writes avgs of each cluster to standard output
def compaction(dfseg,npclusts,outf):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	collection = [] # store compaction calculations
	for x in range(0,len(npclusts)):
		selected = st[npclusts[x]] # constrain seg table to only
		totnps = len(npclusts[x])
		compactclust = (selected.sum(axis=1).values) / totnps
		print(compactclust)
		collection.append(compactclust) # calculate sum, normalized by total # of NPs in the cluster
		out1 = selected.sum(axis=1)
		out1.to_csv("./cluster"+str(x+1)+"-unnormalized.csv")
		out2 = out1 / totnps
		out2.to_csv("./cluster"+str(x+1)+"-normalized.csv")
		print("-- Cluster",str(x+1),"has compaction average",sum(compactclust) / len(compactclust))
	fig, ax = plt.subplots(figsize=(10,10),facecolor='white')
	ax.boxplot(collection,labels=["Cluster " + str(x+1) for x in range(0,len(npclusts))])
	plt.title("Boxplot of Clustered Compactions")
	plt.savefig(outf + "cluster-compaction-boxplot.png")
	print("-- Finished generating boxplot of cluster compactions, saved to",outf + "cluster-compaction-boxplot.png")

def heatmapclust(dfseg,outf):
	st = dfseg.loc[:,(dfseg != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	st = st.T # make NPs be the rows
	fig, ax = plt.subplots(figsize=(10,10),facecolor="white")
	myclust = sns.clustermap(st,metric="jaccard",col_cluster=False)
	plt.gcf().subplots_adjust(bottom=0.20)
	plt.suptitle("Clustered Heatmap of NPs")
	plt.savefig(outf + "cluster-heatmap.png")
	print("-- Finished generating clustered heatmap, saved to",outf + "cluster-heatmap.png")