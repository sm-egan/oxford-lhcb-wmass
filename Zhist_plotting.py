import sys
import numpy as np
import ROOT as R
import csv
import matplotlib.pyplot as plt
from decimal import *
import matplotlib.pyplot as plt
from plot_testing import hist_ratio_plot, simple_scatter
    

if __name__ == "__main__":
    pTmethods = ['GausSmear', 'GausSmear_pTdependent', 'ConstFactor', 'CurveOffset']
    npTmethods = len(pTmethods)

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = './rootfiles/Zsampletest.root'
    input = R.TFile.Open(filename, 'UPDATE')

    if len(sys.argv) > 2:
        ntoys = sys.argv[2]
    else:
        ntoys = 6

    chi2file = './chi2results/Zsampletest_lhcbcuts.csv'  
    firstrow = True
    chi2results = [[]]
    pTParams = [[]]
    with open(chi2file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC)
        for row in readCSV:
            if (firstrow): 
                chi2results[0] = row
                firstrow = False
            else:
                chi2results.append(row)
    
    len(chi2results)
    len(pTParams)

    firstrow = True
    parameterfile = './pTparameters/pTparameters.csv'  
    with open(parameterfile) as csvfile:
        readCSV = csv.reader(csvfile, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC)
        for row in readCSV:
            if (firstrow): 
                pTParams[0] = row
                firstrow = False
            else:
                pTParams.append(row)

    for method in range(0, npTmethods):        
        plot_title = 'Zmumu_' + pTmethods[method] + 'chi2_vs_param'
        simple_scatter(pTParams[method], chi2results[method], plot_title, 'Adjustment parameter', '$\chi^2$ test statistic')
        '''
        for toyit in range (0, ntoys):
            targetH_name = pTmethods[method] + 'Z' + str(toyit)
            hist_ratio_plot(input, 'Znominal', targetH_name)
            plt.close('all')
        '''
