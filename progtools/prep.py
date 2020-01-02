import pandas as pd

# SELECTGRI: takes as input the segmentation table and genomic region of interest file
#	outputs a segmentation table 
def selectgri(stpath,gri):
	fulltable = pd.read_csv(stpath, sep='\t')
	fulltable.insert(loc=0, column="combname", value="chr" + fulltable["chrom"] + ":" + fulltable["start"].map(str) + "-" + fulltable["stop"].map(str))
	st = pd.DataFrame(columns=list(fulltable.columns))
	for ind,row in gri.iterrows(): # get segmentation table rows from each section of the genomic region of interest
		subregion = fulltable[fulltable["chrom"] == "chr13"]
		subregion = subregion[(subregion["start"] >= row[1]) & (subregion["stop"] <= row[2])]
		st = pd.concat([st,subregion]) # add new subregion to result
		st = st.drop_duplicates('combname') # ensure no duplicates
	return st