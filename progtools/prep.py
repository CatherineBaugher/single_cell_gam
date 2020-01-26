'''
PREP.PY
	This file contains initial steps of pipeline,
	including segmentation table processing and generating basic statistics
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# SELECTGRI: takes as input path to segmentation table and genomic region of interest file
#	outputs a segmentation table (with all NPs) constrained to the genomic region of interest
def selectgri(stpath,gri):
	fulltable = pd.read_csv(stpath,sep='\t',low_memory=False)
	fulltable.insert(loc=0,column="window",value=fulltable["chrom"] + ":" + fulltable["start"].map(str) + "-" + fulltable["stop"].map(str))
	st = pd.DataFrame(columns=list(fulltable.columns))
	for ind,row in gri.iterrows(): # get segmentation table rows from each section of the genomic region of interest
		subregion = fulltable[fulltable["chrom"] == "chr13"]
		subregion = subregion[(subregion["start"] >= row[1]) & (subregion["stop"] <= row[2])]
		st = pd.concat([st,subregion]) # add new subregion to result
		st = st.drop_duplicates('window') # ensure no duplicates
	# clean up, so that index is the window name and columns chrom, start and stop are dropped
	st = st.drop(["chrom","start","stop"],axis=1) # keep only the combined window name
	st = st.set_index("window")
	st = st.astype(int)
	return st

# BASICCOUNTS: takes as input the original segmentation table confined to the GRI
#	reports the number of NPs that capture GRI and the number of windows in GRI captured by NPs via stdout
def basiccounts(st):
	# calculate number of NPs that capture GRI
	totnum = len(st.columns)
	onlyrelevant = st.loc[:,(st != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	numcapture = len(onlyrelevant.columns)
	perc = numcapture / totnum * 100
	print("--",str(numcapture),"out of",str(totnum),"NPs capture at least one window in the genomic region of interest","(" + str(round(perc,2)) + "%)")
	# number of windows in GRI captured by NPs
	totwind = len(st.index) # total number of windows
	onlyrelevant = st.loc[(st != 0).any(axis=1)] # drop windows which have no NPs
	numcapture = len(onlyrelevant.index)
	perc = numcapture / totwind * 100
	print("--",str(numcapture),"out of",str(totwind),"windows of the GRI are captured by at least one NP","(" + str(round(perc,2)) + "%)")

# COUNTST: helper function to sum rows or columns in ST and format into a dataframe output
def countst(myst,myaxis,colname):
	if(myaxis == 1):
		ind = myst.index
	else:
		ind = myst.columns
	mycounts = pd.DataFrame(columns=[colname],index=ind)
	mycounts[colname] = myst.sum(axis=myaxis) # sum of nps for each window
	return mycounts

# CHECKCOVERAGE: takes as input the segmentation table confined to the GRI and a path to direct output to
#	saves two csv files which provide counts of data available to NPs and data available to windows in GRI
def checkcoverage(st,outf):
	# make count for each window
	windcounts = countst(st,1,"Number of NPs which capture it")
	windcounts.index.name = "Genomic Window"
	windcounts.to_csv(outf + "GRI-info-count.csv")
	print("-- Finished count of NPs for each window in GRI, saved to",outf + "GRI-info-count.csv")
	# make count for each NP
	npcounts = countst(st,0,"Number of windows it captures")
	npcounts.index.name = "Nuclear Profile"
	npcounts.to_csv(outf + "NP-info-count.csv")
	print("-- Finished count of windows for each NP, saved to",outf + "NP-info-count.csv")

# HISTOGRAM: takes as input the segmentation table confined to the GRI and a path to direct output to
#	saves a histogram of the number of genomic windows in NPs
def histogram(st,outf):
	onlyrelevant = st.loc[:,(st != 0).any(axis=0)] # drop NPS which do not cover anything in GRI
	npcounts = countst(st,0,"Number of genomic windows observed") # count number of windows per NP
	fig, ax = plt.subplots(figsize=(10,10),facecolor='white')
	plt.hist(npcounts.values, bins='auto')
	plt.title("Histogram of the Number of Windows Captured by NPs")
	plt.ylabel("Number of nuclear profiles")
	plt.savefig(outf + "NP-windowcount-histogram.png")
	print("-- Finished generating histogram of window counts, saved to",outf + "NP-windowcount-histogram.png")

# FILTERNPS: takes as input the segmentation table dataframe
#	outputs the segmentation table without NPs which have below some threshold # of windows
def filternps(st,threshold):
	beforelen = len(st.columns)
	npcounts = countst(st,0,"numwindows") # get number of windows to each NP
	selected = npcounts[npcounts["numwindows"] > int(threshold)].index # filter the sums to threshold
	result = st[selected]
	afterlen = len(result.columns)
	print("--",str(afterlen),"out of",str(beforelen),"NPs selected (" + str((afterlen / beforelen) * 100) + ")")
	return result