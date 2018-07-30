import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt


def hist_ratio_plot(rootfileo, nominalHstr, targetHstr, wHstr = "None"):

    nominalH = rootfileo.Get(nominalHstr)
    targetH = rootfileo.Get(targetHstr)
    
    nominalH.Scale(targetH.Integral()/nominalH.Integral())

    nbinsN = nominalH.GetNbinsX()
    nbinsT = targetH.GetNbinsX()
    countsN = []
    countsT = []
    ratioTN = []
    ratioerrTN = []
    #ratioerrNN = []

    #targetH.Scale(nominalH.Integral(30,50)/targetH.Integral(30,50))
    nomEntries = nominalH.GetEntries() - (nominalH.GetBinContent(0) + nominalH.GetBinContent(nbinsN+1)) 
    tarEntries = targetH.GetEntries() - (targetH.GetBinContent(0) + targetH.GetBinContent(nbinsT+1)) 
    
    print(nomEntries)
    print(tarEntries)
    print(targetH.GetEntries())
    print(targetH.GetBinContent(0))
    print(targetH.GetBinContent(nbinsT+1))

    AdjustEntries = 0
    ratiosumT = 0
    entrycountN = 0
    entrycountT = 0

    if not wHstr == "None":
        print("THIRD HISTOGRAM GIVEN, INITIALIZING VARIABLES")
        WpredH = rootfileo.Get(wHstr)
        WpredH.Scale(targetH.Integral(30,50)/WpredH.Integral(30,50))
        nbinsW = WpredH.GetNbinsX()
        countsW = []

        #WpredH.Scale(nominalH.Integral(30,50)/WpredH.Integral(30,50))
        #WpredEntries = WpredH.GetEntries() - (WpredH.GetBinContent(0) + WpredH.GetBinContent(nbinsW+1))
        ratioWN = []
        ratioerrWN = []
        entrycountW = 0
        ratiosumW = 0

    
    #nominalH.Scale(1/nominalH.Integral(30,50))
    #targetH.Scale(1/targetH.Integral(30,50))

    print(nominalH.Integral(30,50))
    print(targetH.Integral(30,50))
    print("LOOPING OVER BINS")

    for bin in range (1, nbinsN + 1):
    
        bincountN = nominalH.GetBinContent(bin)
        binerrorN = nominalH.GetBinError(bin)
        #print(binerrorN)
        binerrorNN = np.sqrt(2)*(binerrorN**2)/bincountN


        countsN.append(bincountN)
        #ratioerrNN.append(binerrorNN)
        entrycountN += bincountN 
        
        bincountT = targetH.GetBinContent(bin)
        binerrorT = targetH.GetBinError(bin)
        #print(binerrorT)
#print('target bin ' + str(bin) + 'content: ' + str(bincountT)
        #bincountT = bincountT*(nomEntries/tarEntries)
    #print('after normalization: ' +  str(bincountnormed))
        countsT.append(bincountT)
        entrycountT += bincountT

        countratioT = bincountT/bincountN
        errorratioT = np.sqrt((binerrorT/bincountN)**2 + (bincountT*binerrorN/(bincountN**2))**2)
        #print(errorratioT)

        ratioTN.append(countratioT)
        ratioerrTN.append(errorratioT)
        ratiosumT += (bincountT - bincountN)/bincountN
    
        if not wHstr == "None":
            bincountW = WpredH.GetBinContent(bin)
            binerrorW = WpredH.GetBinError(bin)
            #print(binerrorW)
            #bincountW = bincountW*(nomEntries/WpredEntries)
            countsW.append(bincountW)
            entrycountW += bincountW

            countratioW = bincountW/bincountN
            errorratioW = np.sqrt((binerrorW/bincountN)**2 + (bincountW*binerrorN/(bincountN**2))**2)
            #print(errorratioW)

            ratioWN.append(countratioW)
            ratioerrWN.append(errorratioW)
            ratiosumW += (bincountW-bincountN)/bincountN
            

#xlims = nominalH.GetXaxis().GetXbins()
#print(type(xlims))
    print('RATIO SUM: ' + str(ratiosumT))
    print(entrycountN)
    print(entrycountT)

# Append one extra item to match the x axis created by np.linspace
    '''
    countsN.append(1)
    countsT.append(1)
    ratioTN.append(1)
    ratioerrTN.append(0)
    ratioerrNN.append(0)
    '''

    ratioTN = np.array(ratioTN)
    ratioerrTN = np.array(ratioerrTN)
    #ratioerrNN = np.array(ratioerrNN)

    histx = np.linspace(30,50,nbinsN, endpoint=True)
    normline = [1 for i in range(nbinsN)]
    normline = np.array(normline)

#plt.show()

    fig1, (axHist,axRatios) = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[5,3]})
    
    axHist.step(histx, countsN, where='post', label='Nominal W mass (80.40 GeV)', c='k')
    
    if not wHstr == "None":
        print("ADDING NEAREST w TEMPLATE TO PLOT")
        '''
        countsW.append(1)
        ratioWN.append(1)
        '''
        ratioWN = np.array(ratioWN)
        #ratioerrWN.append(0)
        ratioerrWN = np.array(ratioerrWN)
        pluserrWN = ratioWN + ratioerrWN
        minuserrWN = ratioWN - ratioerrWN

        axRatios.fill_between(histx, ratioWN+ratioerrWN, ratioWN-ratioerrWN, step='post', alpha=0.5, linestyle='-.', color='m')
        axRatios.step(histx, ratioWN, where='post', c='r')

        template_label = 'Predicted W mass (' + wHstr[2:] + ' GeV)'
        axHist.step(histx, countsW, where='post', label = template_label, c='r')


    axHist.step(histx, countsT, where='post', label=targetHstr, c='b')
    
    axHist.set_xlim(30,50)
    #axHist.set_ylim(0,0.05)
    axHist.set_ylabel('Event Count')

    #axRatios.fill_between(histx, normline+ratioerrNN, normline-ratioerrNN, step='post', alpha=0.5, linestyle='--', color='k')
    axRatios.plot(histx, normline, c='k')
    axRatios.fill_between(histx, ratioTN+ratioerrTN, ratioTN-ratioerrTN, step='post', alpha=0.5, linestyle='--', color='c')
    axRatios.step(histx, ratioTN, where='post', c='b')    

    axRatios.set_xlabel('mu_PT')
    axRatios.set_ylabel('Target/Nominal')

    axHist.legend()
    plot_name = 'plots/hist_ratio_plots/hist_ratios_' + targetHstr
    fig1.savefig(plot_name + '.pdf')
    fig1.savefig(plot_name + '.png')

#for bin in range (:
#   bin = nominalH.GetBinContent(i)

#print(len(bincounts))
#print(bincounts[14]) 
    
if __name__ == "__main__":
    filename = "./rootfiles/build_compare_templatesWm.root"
    input = R.TFile.Open(filename, "UPDATE")

    hist_ratio_plot(input, "Wm80.40", "GausSmearWm3", "Wm80.50")
