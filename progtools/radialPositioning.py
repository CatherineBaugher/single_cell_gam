# These libraries are all needed, and need to be installed if not already installed.
import numpy as np # pip install numpy
import operator
from operator import itemgetter
import statistics as stat # pip install statistics
import seaborn as sns # pip install seaborn
import matplotlib.pyplot as plt # pip install matplotlib
import matplotlib.ticker as ticker
from collections import Counter
import pandas as pd
# This is only used to get rid of the warnings that do not necessarily mean anything, you can comment the next 2 lines if you want.
# import warnings
# warnings.filterwarnings("ignore")


# This class is all of the functions that I used to calculate radial positioning within a cell.
class radialPosition:
    # This is the initializer. This function takes all of the parameters so that they only need to be passed once.
    # NPNames is a list of all Nuclear Profile Names
    # WindowSums will be a summation of all data within a row, or sumsrow from my previous implementation
    # NPSums will be a summation of all data within a Column, or sumscol from previous implementation
    def __init__(self):
        pass

    # radialPositionCalc takes only the apical and equatorial cutoff values as an integer. 
    # For Example, if you want an apical cutoff of 25% and an equatorial cutoff of 75%, you
    # need to run the function as follows: radialPositionCalc(25, 75)
    #
    # This function will return a list with the following items in the given order:
    #   FinalAbove = list of Nuclear Profiles and the number windows they appear in ascending order for the Equatorial Data
    #   FinalABelow = list of Nuclear Profiles and the number windows they appear in ascending order for the Apical Data
    #   Above = The cutoff Window amount of the Equatorial Data
    #   Below = The cutoff Window amount of the Apical Data
    #   len(FinalAboveIndices) = number of Equatorial Nuclear Profiles
    #   len(FinalBelowIndices) = number of Apical Nuclear Profiles
    
    def radialPositionCalc(self, inputdata, apical=25, equatorial=75):

        inputdata = inputdata[inputdata["Number of windows it captures"] > 0]

        NPSums = inputdata["Number of windows it captures"]
        NPNames = inputdata.index

        # print(NPSums)
        # print(NPNames)

        if inputdata.empty:
            return []


        # Calculating the cutoff for the given quartiles.
        Below = np.percentile(NPSums, apical)
        Above = np.percentile(NPSums, equatorial)

        # Getting the Indices of the elements that are above and below the quartiles. This is used to 
        # generate the set of features that are either Apical or Equatorial.
        BelowIndices = np.argwhere(NPSums <= Below)
        AboveIndices = np.argwhere(NPSums >= Above)

        FinalBelowIndices = []
        FinalAboveIndices = []

        # Fixing a formatting issue.
        for i in range(0, BelowIndices.shape[0]):
            FinalBelowIndices.append(BelowIndices[i][0])

        # Fixing a formatting issue
        for i in range(0, AboveIndices.shape[0]):
            FinalAboveIndices.append(AboveIndices[i][0])


        BelowHelper = []
        AboveHelper = []

        # Generating a helper list so that I can sort the data later. Helper list is made up of 
        # elements in the following format: (Feature Name, Number of Windows that feature appears in.)
        for i in range(0, len(FinalBelowIndices)):
            BelowHelper.append((NPNames[FinalBelowIndices[i]], NPSums[FinalBelowIndices[i]]))

        # Sorting the list of elements for more friendly view.
        FinalBelow = sorted(BelowHelper, key=itemgetter(1), reverse=False)

        # Generating a helper list so that I can sort the data later. Helper list is made up of 
        # elements in the following format: (Feature Name, Number of Windows that feature appears in.)
        for i in range(0, len(FinalAboveIndices)):
            AboveHelper.append((NPNames[FinalAboveIndices[i]], NPSums[FinalAboveIndices[i]]))

        # Sorting the list of elements for more friendly view.
        FinalAbove = sorted(AboveHelper, key=itemgetter(1), reverse=False)

        return(FinalAbove, FinalBelow, Above, Below, len(FinalAboveIndices), len(FinalBelowIndices))


    # clusterBoxPlot
    # This function generates the box plot based on the nuclear profile information.
    #
    # Input: the output directory (2 files are loaded from the output directory, NP-info-count.csv and cluster-kmeans-raw.csv are both
    # required to run the program, and if the files do not exist, there will be errors.)
    #
    # Output: there will be a new file saved called NP-ClusterBoxPlot.png within the output directory that is a boxplot of the cluster
    # data. There is nothing returned and nothing printed to standard output currently.
    def clusterBoxPlot(self, outdir):
        NPInfo = outdir + "NP-info-count.csv"
        clusterInfo = outdir + "cluster-kmeans-raw.csv"

        data = pd.read_csv(clusterInfo, header=None)
        alldata = pd.read_csv(NPInfo)

        allNames = alldata.iloc[:,0].values
        WindowCount = alldata.iloc[:,1].values

        NPNames = data.iloc[1:,0].values
        ClusterIndex = data.iloc[1:,1].values

        cluster1Indices = []
        cluster2Indices = []
        cluster3Indices = []

        cluster1Indices = np.argwhere(ClusterIndex == '0')
        cluster2Indices = np.argwhere(ClusterIndex == '1')
        cluster3Indices = np.argwhere(ClusterIndex == '2')

        cluster1Names = NPNames[cluster1Indices]
        cluster2Names = NPNames[cluster2Indices]
        cluster3Names = NPNames[cluster3Indices]

        cluster1FinalIndices = []
        cluster2FinalIndices = []
        cluster3FinalIndices = []

        for i in cluster1Names:
            cluster1FinalIndices.append(np.argwhere(allNames == i)[0][:][0])

        for i in cluster2Names:
            cluster2FinalIndices.append(np.argwhere(allNames == i)[0][:][0])

        for i in cluster3Names:
            cluster3FinalIndices.append(np.argwhere(allNames == i)[0][:][0])

        cluster1FinalValues = WindowCount[cluster1FinalIndices]
        cluster2FinalValues = WindowCount[cluster2FinalIndices]
        cluster3FinalValues = WindowCount[cluster3FinalIndices]

        data_to_plot = [cluster1FinalValues, cluster2FinalValues, cluster3FinalValues]

        fig = plt.figure(1)
        plt.xlabel('Clusters')
        plt.ylabel('Windows')

        ax = fig.add_subplot()
        ax.set_xticklabels(['Cluster1', 'Cluster2', 'Cluster3'])

        bp = ax.boxplot(data_to_plot, patch_artist=True)

        for box in bp['boxes']:
            box.set( color='#7570b3', linewidth=2)
            box.set( facecolor = '#1b9e77' )

        for whisker in bp['whiskers']:
            whisker.set(color='#7570b3', linewidth=2)

        for cap in bp['caps']:
            cap.set(color='#7570b3', linewidth=2)

        for median in bp['medians']:
            median.set(color='#b2df8a', linewidth=2)

        for flier in bp['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.5)

        fig.savefig((outdir + "NP-ClusterBoxPlot.png"), bbox_inches='tight', pad_inches=0.5)

        return