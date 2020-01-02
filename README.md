# Command Line Toolkit for Single-Cell GAM Analysis
#### Table of Contents
[Setup](#setup)  
[Obtaining Data](#data)
[Running](#run)

<a name="setup"/>

## Setup
### Python 3
Install the latest version of Python 3.

[Anaconda](https://www.anaconda.com/distribution/) may also be installed for easier setup of required packages.

### Required Packages
	* [matplotlib](https://matplotlib.org/users/installing.html)
	* [pandas](https://pandas.pydata.org/pandas-docs/stable/install.html)

### Installation
Save the github repository via `git clone https://github.com/CatherineBaugher/single_cell_gam.git`.

To ensure a proper Python environment, a test should be performed using the command `python3 singlecellgam.py ./exampledata/chr13seg.txt ./exampledata/hist1.bed -h`.

This should display the --help information:
```
usage: singlecellgam.py [-h] [-b] [-s] [-p] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  -b, --basicstats      Performs the following basic statistics, using as
                        input segregation table and GRI file: -- Count the
                        number of NPs which capture the GRI -- Count the
                        number of windows in the GRI which are captured by NPs
                        in the segregation table -- For each window within the
                        GRI, calculate the number of NPs that detected it
  -s, --similarity      Generate a normalized similarity matrix from a given
                        segregation table and GRI file
  -p, --pca             Generate a PCA from a given segregation table and GRI
                        file
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        Specify a directory to save any outputted files to
```

<a name="data"/>

## Obtaining Data
### Segmentation Table
The full 30kb GAM segmentation table may be downloaded directly at: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE64881&format=file&file=GSE64881%5Fsegmentation%5Fat%5F30000bp%2Epassqc%2Emultibam%2Etxt%2Egz

This file and a 1 Mb resolution segmentation table are also hosted on the GEO archive: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881

### Genomic Region of Interest (GRI)
A file indicating the region of interest to analyze is required. This BED file should have 3 fields:
	* Chromosome (chr#)
	* Start position of the region
	* Stop position of the region
The file may contain multiple lines. This may allow the user to, for example, skip unmappable regions.

<a name="run"/>

## Running
### Arguments
[placeholder]

### Flags
[placeholder]

### Example Tutorial
[placeholder]