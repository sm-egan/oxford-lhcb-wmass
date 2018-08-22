import os 
import sys
import numpy as np
import ROOT as R
from plot_testing import hist_ratio_plot

if __name__ == "__main__":
    if len(sys.argv) > 1:
        rootfile = sys.argv[1]
    else:
        rootfile = 'rootfiles/Upsilonsample_20082018-1440.root'
    input = R.TFile.Open(rootfile, "UPDATE")


    if len(sys.argv) == 3:
        nominalH = sys.argv[2]
        targetH = sys.argv[3]
        if nominalH.find('deltamuPT') > -1:
            xmin = -20
            xmax = 20
        elif nominalH.find('ratiomuPT') > -1:
            xmin = 0
            xmax = 10            
        hist_ratio_plot(input, nominalH, targetH, xmin, xmax)

    else:
        curveobs = ['deltamuPT', 'deltamuP', 'asym_dmuPT', 'asym_dmuP']
        
        for coindex in range(0,len(curveobs)):
            print('Preparing histograms for ' + curveobs[coindex])
            nominalH = curveobs[coindex] + '_CurveOffsetUpsilon0'

            if nominalH.find('asym') > -1:
                xmin = -1
                xmax = 1
            elif nominalH.find('deltamuP') > -1:
                xmin = -20
                xmax = 20
            elif nominalH.find('ratiomuPT') > -1:
                xmin = 0
                xmax = 10            

            for toyindex in range(1,6):
                print('TOY ' + str(toyindex))
                targetH = curveobs[coindex] + '_CurveOffsetUpsilon' + str(toyindex)
                hist_ratio_plot(input, nominalH, targetH, xmin, xmax)
                print('plot saved')
