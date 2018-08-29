import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt
import sys

def simple_scatter(x, y, title='', xlabel='x', ylabel='y'):
    fig, ax = plt.subplots(1);
    ax.plot(x, y, 'bo')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if title == '':
        title = 'untitledplot'

    fig.savefig('plots/' + title + '.png')
    fig.savefig('plots/' + title + '.pdf')

def 2d_chi2_pTparams ():


def plot_root_hist (rootfileo, Hstr):
    H = rootfileo.Get(Hstr)
    nbins = H.GetNbinsX()
    counts = []
    errors = []

    for bin in range (1, nbins+1):
        bincount = H.GetBinContent(bin)
        binerror = H.GetBinError(bin)

        counts.append(bincount)
        errors.append(binerror)

    fig, ax = plt.subplots(1)
    if Hstr.find('Z') > -1:
        label = 'Nominal Z'
        xmin = 80
        xmax = 100
        xlabel = 'Dimuon invariant mass (GeV)'
    elif Hstr.find('Upsilon') > -1:
        nominal_label = 'Nominal Upsilon'
        xmin = 9.45601
        xmax = 9.45605
        xlabel = 'Dimuon invariant mass (GeV)'
    else:
        nominal_label = 'Nominal W mass (80.40 GeV)'
        xmin = 30
        xmax = 50
        xlabel = 'Muon pT (GeV)'
    

    histx = np.linspace(xmin, xmax, nbins, endpoint=True)
    ax.step(histx, counts, where='post', c='k')   
    ax.errorbar(histx, counts, yerr=errors, fmt='none',capsize=3)
    
    ax.set_xlim(xmin, xmax)
    #axHist.set_ylim(0,0.05)
    ax.set_ylabel('Event Count')

    #axRatios.fill_between(histx, normline+ratioerrNN, normline-ratioerrNN, step='post', alpha=0.5, linestyle='--', color='k')
    ax.set_xlabel(xlabel)
    ax.set_title(Hstr)
    
    '''
    if nominalHstr.find('Z') > -1: 
        axHist.legend(loc = 'upper center', bbox_to_anchor=(0., 1.02, 1., .102), ncol=legendcol, mode = 'expand', borderaxespad=0.)
    else:
        axHist.legend()
    '''
    plot_name = 'plots/pyplot_histograms/' + Hstr
    fig.savefig(plot_name + '.pdf', bbox_inches='tight')
    fig.savefig(plot_name + '.png', bbox_inches='tight')
    

def hist_ratio_plot(rootfileo, nominalHstr, targetHstr, wHstr = 'None', xmin = "", xmax = ""):

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

    if wHstr.find('None') ==  -1:
        print("THIRD HISTOGRAM GIVEN, INITIALIZING VARIABLES")
        WpredH = rootfileo.Get(wHstr)
        WpredH.Scale(targetH.Integral()/WpredH.Integral())
        nbinsW = WpredH.GetNbinsX()
        countsW = []

        #WpredH.Scale(nominalH.Integral(30,50)/WpredH.Integral(30,50))
        #WpredEntries = WpredH.GetEntries() - (WpredH.GetBinContent(0) + WpredH.GetBinContent(nbinsW+1))
        ratioWN = []
        ratioerrWN = []
        entrycountW = 0
        ratiosumW = 0

    
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
        
        if bincountN == 0 or bincountT == 0:
            countratioT = 1
            errorratioT = 0            
        else:
            countratioT = bincountT/bincountN
            errorratioT = np.sqrt((binerrorT/bincountN)**2 + (bincountT*binerrorN/(bincountN**2))**2)
        
        #print(errorratioT)

        ratioTN.append(countratioT)
        ratioerrTN.append(errorratioT)
        #ratiosumT += (bincountT - bincountN)/bincountN
    
        if wHstr.find('None') == -1:
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

    #xmin = nominalH.GetXaxis().GetXmin()
    #xmax = nominalH.GetXaxis().GetXmin()
    
    normline = [1 for i in range(nbinsN)]
    normline = np.array(normline)

#plt.show()

    fig1, (axHist,axRatios) = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[5,3]})
    
    #You should find a way to get these numbers automatically
    if nominalHstr.find('Z') > -1:
        nominal_label = 'Nominal Z'
        if xmin == "" or xmax == "":
            xmin = 80
            xmax = 100
        xlabel = 'Dimuon invariant mass (GeV)'
    elif nominalHstr.find('Upsilon') > -1:
        nominal_label = 'Nominal Upsilon'
        if xmin == "" or xmax == "":
            xmin = 9.45601
            xmax = 9.45605
        if nominalHstr.find('deltamuPT') > -1:
            xlabel = '$p_T^{\mu^-} - p_T^{\mu^+}$ ($GeV/c$)'
        elif nominalHstr.find('deltamuP') > -1:
            xlabel = '$p^{\mu^-} - p^{\mu^+}$ ($GeV/c$)'
        elif nominalHstr.find('asym_dmuP') > -1:
            xlabel = '$p^{\mu^-} - p^{\mu^+}$ / $p^{\mu^-} + p^{\mu^+}$'
        elif nominalHstr.find('asym_dmuPT') > -1:
            xlabel = '$p_T^{\mu^-} - p_T^{\mu^+}$ / $p_T^{\mu^-} + p_T^{\mu^+}$'
    else:
        nominal_label = 'Nominal W mass (80.40 GeV)'
        if xmin == "" or xmax == "":
            xmin = 30
            xmax = 50
        xlabel = 'Muon pT (GeV)'

    histx = np.linspace(xmin, xmax, nbinsN, endpoint=True)
    axHist.step(histx, countsN, where='post', label=nominal_label, c='k')   
    legendcol = 2
    
    if wHstr.find('None') == -1:
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
        
        axRatios.errorbar(histx, ratioWN, yerr=ratioerrWN, fmt='none',capsize=3, color = 'm')
        #axRatios.fill_between(histx, ratioWN+ratioerrWN, ratioWN-ratioerrWN, step='post', alpha=0.5, linestyle='-.', color='m')
        axRatios.step(histx, ratioWN, where='post', c='r')

        pred_mass = (int(wHstr[-1]) - 6)*0.1 + 80.40
        template_label = 'Predicted W mass (' +  str(pred_mass) + ' GeV)'
        axHist.step(histx, countsW, where='post', label = template_label, c='r')
        legendcol = 3

    axHist.step(histx, countsT, where='post', label=targetHstr, c='b')
    
    axHist.set_xlim(xmin, xmax)
    #axHist.set_ylim(0,0.05)
    axHist.set_ylabel('Event Count')

    #axRatios.fill_between(histx, normline+ratioerrNN, normline-ratioerrNN, step='post', alpha=0.5, linestyle='--', color='k')
    axRatios.plot(histx, normline, c='k')
    axRatios.errorbar(histx, ratioTN, yerr=ratioerrTN, fmt='none',capsize=3, color='c')
    #axRatios.fill_between(histx, ratioTN+ratioerrTN, ratioTN-ratioerrTN, step='post', alpha=0.5, linestyle='--', color='c')
    axRatios.step(histx, ratioTN, where='post', c='b')    

    axRatios.set_xlabel(xlabel)
    axRatios.set_ylabel('Toy/Template')
    
    if nominalHstr.find('Z') or nominalHstr.find('Upsilon') > -1: 
        axHist.legend(loc = 'upper center', bbox_to_anchor=(0., 1.02, 1., .102), ncol=legendcol, mode = 'expand', borderaxespad=0.)
    else:
        axHist.legend()

    plot_name = 'plots/hist_ratio_plots/hist_ratios_' + targetHstr
    fig1.savefig(plot_name + '.pdf', bbox_inches='tight')
    fig1.savefig(plot_name + '.png', bbox_inches='tight')

#for bin in range (:
#   bin = nominalH.GetBinContent(i)

#print(len(bincounts))
#print(bincounts[14]) 
    
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "./rootfiles/Zsampletest.root"
    input = R.TFile.Open(filename, "UPDATE")
    
    if filename.find('Z') > -1:
        nominalH = 'Znominal'
        prop = 'Z'
    elif filename.find('Wm') > -1:
        nominalH = 'Wm80.40'
        prop = 'Wm'
    else:
        nominalH = 'Wp80.40'
        prop = 'Wp'

    targetH = 'GausSmear' + prop + '0'

    if len(sys.argv) == 3:
        H = sys.argv[2]
        plot_root_hist(input, H)

    elif len(sys.argv) == 4:
        nominalH = sys.argv[2]
        targetH = sys.argv[3]
        hist_ratio_plot(input, nominalH, targetH)

    elif len(sys.argv) > 4:
        nominalH = sys.argv[2]
        targetH = sys.argv[3]
        compareH = sys.argv[4]
        hist_ratio_plot(input, nominalH, targetH, compareH)
    
