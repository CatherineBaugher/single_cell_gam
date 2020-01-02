import argparse
import pandas as pd
import progtools.prep
import progtools.variation

parser = argparse.ArgumentParser()
# ARGUMENTS
parser.add_argument("st", help="Input tab-separated GAM segregation table")
parser.add_argument("gri", help="Input bed file of genomic region of interest, with each line of the form chr# | start | end")
# FLAGS
parser.add_argument("-b","--basicstats", help="Performs the following basic statistics, using as input segregation table and GRI file:"\
		"\n-- Count the number of NPs which capture the GRI" \
		"\n-- Count the number of windows in the GRI which are captured by NPs in the segregation table" \
		"\n-- For each window within the GRI, calculate the number of NPs that detected it",action="store_true")
parser.add_argument("-s","--similarity", help="Generate a normalized similarity matrix from a given segregation table and GRI file",action="store_true")
parser.add_argument("-p","--pca", help="Generate a PCA from a given segregation table and GRI file",action="store_true")
# OPTIONS
parser.add_argument("-o","--outputfile", help="Specify a directory to save any outputted files to")

# Parse input and call necessary processing functions
args = parser.parse_args()
dfgri = pd.read_csv(args.gri, header=None, sep='\t')
segtab = progtools.prep.selectgri(args.st,dfgri) # Specifically select data from segmentation table within the genomic region of interest
if args.basicstats:
	print("basic stats chosen")
if args.similarity:
	print("similarity chosen")
if args.pca:
	print("pca chosen")