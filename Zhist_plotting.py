import sys
import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt
from plot_testing import hist_ratio_plot
    

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

    for method in range(0, npTmethods):
        
        for toyit in range (0, ntoys):
            targetH_name = pTmethods[method] + 'Z' + str(toyit)
            hist_ratio_plot(input, 'Znominal', targetH_name)
            plt.close('all')
