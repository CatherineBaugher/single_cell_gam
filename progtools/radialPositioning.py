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
import warnings
warnings.filterwarnings("ignore")


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