import sys
import numpy as np
import ROOT as R
import csv
import matplotlib.pyplot as plt
from decimal import *
import matplotlib.pyplot as plt
from plot_testing import hist_ratio_plot, simple_scatter
    

if __name__ == "__main__":
    pTmethods = ['GausSmear', 'GausSmear_pTdependent', 'MomentumScale', 'CurvatureBias']
    npTmethods = len(pTmethods)

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = './rootfiles/Zsampletest.root'
    input = R.TFile.Open(filename, 'UPDATE')

    ntoys = 11

    if len(sys.argv) > 2:
        chi2file = sys.argv[2]
    else:
        chi2file = './chi2results/Zsampletest_lhcbcuts.csv'  

    firstrow = True
    chi2results = [[]]
    pTParams = [[]]
    chi2count = 0
    with open(chi2file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC)
        for row in readCSV:
            if chi2count > 3:
                break 
            else:
                if (firstrow): 
                    chi2results[0] = row
                    firstrow = False
                else:
                    chi2results.append(row)
                chi2count += 1
    
    len(chi2results)
    len(pTParams)

    firstrow = True
    if filename.find('Upsilon') > -1:
        propagator = 'Upsilon'
        nominalH_name = 'Upsilontemplate'
    else:
        propagator = 'Z'
        nominalH_name = propagator + 'nominal'

    if len(sys.argv) > 3:
        parameterfile = sys.argv[3]
    else:
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
        plot_title = propagator + 'mumu_' + pTmethods[method] + 'chi2_vs_param'
        simple_scatter(pTParams[method], chi2results[method], plot_title, 'Adjustment parameter', '$\chi^2$ test statistic')
        
        for toyit in range (0, ntoys):
            targetH_name = pTmethods[method] + propagator + str(toyit)
            print('Attempting hist_ratio plot of ' + targetH_name + ' and ' +nominalH_name)
            hist_ratio_plot(input, nominalH_name, targetH_name)
            plt.close('all')
        
