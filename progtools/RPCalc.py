import progtools.radialPositioning
import progtools.prep
import pandas as pd

def RPCall(dfseg, outdir, Apical=25, Equatorial=75):
    # tmp = []
    rp = progtools.radialPositioning.radialPosition()
    if progtools.prep.checkcoverage(dfseg, outdir) == None:
        print("-- No Data to Calculate.")
        return
    tmp = rp.radialPositionCalc(progtools.prep.checkcoverage(dfseg, outdir)[1], Apical, Equatorial)


    if tmp == []:
        print("-- No Data to Calculate.")
        return

    ENames = []
    EData = []
    ANames = []
    AData = []

    for i in range(0,len(tmp[0])):
        ENames.append(tmp[0][i][0])
        EData.append(tmp[0][i][1])

    for i in range(0,len(tmp[1])):
        ANames.append(tmp[1][i][0])
        AData.append(tmp[1][i][1])



    EquatorialDF = pd.DataFrame(EData, index=ENames, columns = ['Window Count'])
    ApicalDF = pd.DataFrame(AData, index=ANames, columns = ['Window Count'])

    fileLocE = outdir + 'Equatorial.csv'
    fileLocA = outdir + 'Apical.csv'

    EquatorialDF.to_csv(fileLocE)
    ApicalDF.to_csv(fileLocA)

    print("-- Number of Equatorial Elements:", tmp[4])
    print("-- Number of Apical Elements:", tmp[5])
    print("-- Equatorial Cutoff:", int(tmp[2]))
    print("-- Apical Cutoff:", int(tmp[3]))