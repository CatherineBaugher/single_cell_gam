# Command Line Toolkit for Single-Cell GAM Analysis
#### Table of Contents
* [Setup](#setup)
  * [Python 3](#python)
  * [Required Packages](#packages)
  * [Installation](#install)
* [Obtaining Data](#data)
  * [Segmentation Table](#st)
  * [GRI File](#gri)
* [Running](#run)
  * [Arguments](#arguments)
  * [Flags](#flags)
* [Example Tutorial](#example)

<br><br><br>

<a name="setup"/>

## Setup

<a name="python"/>

### Python 3
Install the latest version of [Python 3](https://www.python.org/downloads/).

[Anaconda](https://www.anaconda.com/distribution/) may also be installed for easier setup of required packages.

<a name="packages"/>

### Required Packages
* [matplotlib](https://matplotlib.org/users/installing.html)
* [pandas](https://pandas.pydata.org/pandas-docs/stable/install.html)

<a name="install"/>

### Installation
No installation is required. Simply save the github repository via `git clone https://github.com/CatherineBaugher/single_cell_gam.git`.

To ensure a proper Python environment, perform a test using the command `python3 singlecellgam.py -h` within the root directory of the project.

This should display the --help information:
```
usage: singlecellgam.py [-h] [-b] [-v] [-c] [-o OUTPUTDIR] st gri

positional arguments:
  st                    Input tab-separated GAM segregation table
  gri                   Input bed file of genomic region of interest, with
                        each line of the form chr# | start | end

optional arguments:
  -h, --help            show this help message and exit
  -b, --basicstats      Performs the following basic statistics, using as
                        input segregation table and GRI file: -- Count the
                        number of NPs which capture the GRI -- Count the
                        number of windows in the GRI which are captured by NPs
                        in the segregation table -- For each window within the
                        GRI, calculate the number of NPs that detected it
  -v, --variation       Performs the following methods to assess variation: --
                        Generate a normalized similarity matrix from a given
                        segregation table and GRI file -- Generate a PCA from
                        a given segregation table and GRI file
  -c, --cluster         Performs the following clustering methods: -- Cluster
                        heatmap of NPs that captured similar subsets of
                        windows from the GRI -- Cluster correlation matrix of
                        NPs
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Specify a directory to save any outputted files to
```

<a name="data"/>

## Obtaining Data

<a name="st"/>

### Segmentation Table
The full 30kb GAM segmentation table may be downloaded directly at: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE64881&format=file&file=GSE64881%5Fsegmentation%5Fat%5F30000bp%2Epassqc%2Emultibam%2Etxt%2Egz

This file and a 1 Mb resolution segmentation table are also hosted on the [GEO archive](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881).

<a name="gri"/>

### Genomic Region of Interest (GRI)
A file indicating the region of interest to analyze is required. This BED file should have 3 fields:
1. Chromosome (chr#)
2. Start position of the region
3. Stop position of the region

The GRI file may contain multiple lines. This allows the user to, for example, skip unmappable regions.

<a name="run"/>

## Running

<a name="arguments"/>

### Arguments
There are two required arguments which are required in order to run the script. In order:
* **st**: Path to segmentation table
* **gri**: Path to GRI file

Additionally, an optional argument may be passed with **--outputdir**, indicating a path to create a folder to direct output to. Any files that the program generates will be saved there. By default, files will simply output within the `single_cell_gam` directory.

<a name="flags"/>

### Flags
Flags are used to indicate which parts of the pipeline you wish to run.
* **-b**, **--basicstats**: An initial exploratory analysis is performed to insure that the GRI is adequately captured by NPs in the GAM segregation table. Some basic statistics are recorded to standard output and two files indicating the amount of information within GRI windows and available NPs are generated.
* **-v**, **--variation**: In order to assess variation among NPs relevant to the GRI, a similarity matrix (heatmap) is generated, as well as a 2D linear PCA.
* **-c**, **--cluster**: Unique patterns of contacts within the region of interest may be discovered through clustering. A heatmap of NPs that captured similar subsets of windows from the GRI and a correlation matrix of NPs is generated. These matrices are then clustered by calculating a distance metric, defined as "the proportion of bits in which only one is on amongst those in which at least one is on" [[binary dist as defined by R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist)].

<a name="example"/>

### Example Tutorial
[placeholder]

### TODO
[✓] Create skeleton of argparse program

[✓] Implement basic statistics

[_] Implement variation analysis

[_] Implement cluster generation

[_] Investigate relationship of R 'binary' distance to Python 'jaccard' distance

[_] Write example tutorial
