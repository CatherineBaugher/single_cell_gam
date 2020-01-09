'''
PREP.PY
	This file contains initial steps of pipeline,
	including segmentation table processing and generating basic statistics
'''

import pandas as pd

# SELECTGRI: takes as input path to segmentation table and genomic region of interest file
#	outputs a segmentation table (with all NPs) constrained to the genomic region of interest
def selectgri(stpath,gri):
	fulltable = pd.read_csv(stpath,sep='\t')
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
	perc = numcapture / totnum * 100
	print("--",str(numcapture),"out of",str(totnum),"windows of the GRI are captured by at least one NP","(" + str(round(perc,2)) + "%)")

# CHECKCOVERAGE: takes as input the segmentation table confined to the GRI and a path to direct output to
#	saves two csv files which provide counts of data available to NPs and data available to windows in GRI
def checkcoverage(st,outf):
	# make count for each window
	windcounts = pd.DataFrame(columns=["Number of NPs which capture it"],index=st.index)
	windcounts["Number of NPs which capture it"] = st.sum(axis=1) # sum of nps for each window
	windcounts.index.name = "Genomic Window"
	windcounts.to_csv(outf + "GRI-info-count.csv")
	print("-- Finished count of NPs for each window in GRI, saved to",outf + "GRI-info-count.csv")
	# make count for each NP
	npcounts = pd.DataFrame(columns=["Number of windows it captures"],index=st.columns)
	npcounts.index.name = "Nuclear Profile"
	npcounts["Number of windows it captures"] = st.sum(axis=0) # sum of windows for each np
	npcounts.to_csv(outf + "NP-info-count.csv")
	print("-- Finished count of windows for each NP, saved to",outf + "NP-info-count.csv")
	return (windcounts, npcounts)