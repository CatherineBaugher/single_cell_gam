import argparse
import pandas as pd
from pathlib import Path
import progtools.prep
import progtools.variation
import progtools.cluster

parser = argparse.ArgumentParser()
# ARGUMENTS
parser.add_argument("st", help="Input tab-separated GAM segregation table")
parser.add_argument("gri", help="Input bed file of genomic region of interest, with each line of the form chr# | start | end")
# FLAGS
parser.add_argument("-b","--basicstats", help="Performs the following basic statistics, using as input segregation table and GRI file:"\
		"\n-- Count the number of NPs which capture the GRI" \
		"\n-- Count the number of windows in the GRI which are captured by NPs in the segregation table" \
		"\n-- For each window within the GRI, calculate the number of NPs that detected it",action="store_true")
parser.add_argument("-v","--variation", help="Performs the following methods to assess variation:"\
		"\n-- Generate a normalized similarity matrix from a given segregation table and GRI file"\
		"\n-- Generate a PCA from a given segregation table and GRI file",action="store_true")
parser.add_argument("-c","--cluster", help="Performs the following clustering methods:"\
		"\n-- Cluster heatmap of NPs that captured similar subsets of windows from the GRI"\
		"\n-- Cluster correlation matrix of NPs",action="store_true")
# OPTIONS
parser.add_argument("-o","--outputdir", help="Specify a directory to save any outputted files to")

# parse input and call necessary processing functions
args = parser.parse_args()
dfgri = pd.read_csv(args.gri, header=None, sep='\t')
dfseg = progtools.prep.selectgri(args.st,dfgri) # specifically select data from segmentation table within the genomic region of interest
outdir = "./" # by default, output files to the current directory
print("-------------------------------------------")
if args.outputdir != None:
	Path(args.outputdir).mkdir(parents=True, exist_ok=True) # create the directory to output to
	outdir = args.outputdir
	if(outdir[len(outdir)-1] != '/'):
		outdir = outdir + '/' # ensure directory ends in forward slash
if args.basicstats:
	print("Performing BASIC STATISTICS...")
	progtools.prep.basiccounts(dfseg)
	progtools.prep.checkcoverage(dfseg,outdir)
	print("BASIC STATISTICS done!")
	print("-------------------------------------------")
if args.variation:
	print("Performing NP VARIATION ANALYSIS...")
	progtools.variation.similarity(dfseg,outdir)
	progtools.variation.pca(dfseg,outdir)
	print("NP VARIATION ANALYSIS done!")
	print("-------------------------------------------")
if args.cluster:
	print("Performing NP VARIATION ANALYSIS...")
	# clustering functions here
	print("NP VARIATION ANALYSIS done!")
	print("-------------------------------------------")