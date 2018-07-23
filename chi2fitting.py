#from ROOT import gROOT as R
#R.SetBatch(True) 
#import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt

input = R.TFile.Open("./rootfiles/build_compare_templates.root", "UPDATE")

pTParams = [[]]
Wmass_pred = [[]*2]
Wmass_err = [[]*2]
chi2min_fit = [[]*2]

quadfit = R.TF1("chi2_quadfit", "[0]+(1/([2]**2))*(x-[1])**2", 79.8, 80.8)
quadfit.SetParName(0,"MinChi2")
quadfit.SetParName(1,"WMass")
quadfit.SetParName(2,"WMassUncertainty")


firstrow = True 
with open('./pTparameters.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC)
    for row in readCSV:
        if (firstrow): 
            pTParams[0] = row
            firstrow = False
        else:
            pTParams.append(row)

print(pTParams)
print len(pTParams)
print(pTParams[1])


for toyset in range (0,2):
    print('Now plotting set of toys: ' + str(toyset))
    ntoys = len(pTParams[toyset])
    for toyit in range (0, ntoys):    

        chi2Plot_name = "chi2plot" + str(toyset) + str(toyit)    
        print("Loading graph: " + chi2Plot_name)
        chi2Plot = input.Get(chi2Plot_name)
        print(type(chi2Plot))
    
        chi2stats = chi2Plot.GetY()
    
        p0_init = min(chi2stats)
        print("p0_init: " + str(p0_init))
        p1_init = 80.5
        p2_init = 0.01
    
    #print(chi2max)
    #print("Spread of chi2 results: " + str(chi2max-p0_init))
    #p2_init = (chi2max-p0_init)/(79.8-p1_init)
    #print("Initiating parameter 2 as: " + str(p2_init))
    #p2_init = 1/(abs(p2_init))

        quadfit.SetParameters(p0_init, p1_init, p2_init)

        chi2Plot.Fit("chi2_quadfit")

    #chi2min = quadfit.GetMinimum()
        Mw = quadfit.GetParameter(1)
        MwSigma = quadfit.GetParameter(2)

        print('The best W mass is ' + str(Mw) + '+/-' + str(MwSigma))
        
        try:
            Wmass_pred[toyset].append(Mw)
            Wmass_err[toyset].append(MwSigma)
            chi2min_fit[toyset].append(quadfit.GetParameter(0))
        except:
            Wmass_pred.append([Mw])
            Wmass_err.append([MwSigma])
            chi2min_fit.append([quadfit.GetParameter(0)])
            
        
    #pTParams.append(toyit*pTincrement)

        c = R.TCanvas()
        chi2Plot.Draw("AP")

    #TLegend *legend = new TLegend()
    #legend->SetHeader()
    #legend->Draw()

        chi2image = "./plots/chi2Plot" + str(toyset) +  str(toyit) + "_fit.png"
        c.Print(chi2image)
        c.Close()

    '''
    if not ((len(Wmass_pred) == len(Wmass_err)) and (len(Wmass_pred) == len(chi2min_fit)) and (len(Wmass_pred) == len(pTParams))):
        print("Warning!  The arrays you wish to plot are not of equal length!")
        print('Wmass predictions: ' + str(len(Wmass_pred)))
        print('Wmass errors: ' + str(len(Wmass_err)))
        print('pT parameters: ' + str(len(pTParams)))
        print('chi2 test statistics: ' + str(len(chi2min_fit)))
    '''

fig1, (axW, axChi) = plt.subplots(2, sharex=True)
axW.scatter(pTParams[0], Wmass_pred[0])
axW.errorbar(pTParams[0], Wmass_pred[0], yerr=Wmass_err[0], fmt='o')
axW.set_ylabel('Predicted mass of W (GeV)')

max1 = max(pTParams[0]) + 0.1*max(pTParams[0])
min1 = -0.1*max(pTParams[0])
axW.set_xlim([min1,max1])

axChi.scatter(pTParams[0], chi2min_fit[0])
axChi.set_ylabel('Test stat. of W template fit')
axChi.set_xlabel('Mean percent pT smear (Gaussian)')

fig1.savefig("./plots/gaus_pT_smear.png")


fig2, (ax2W, ax2Chi) = plt.subplots(2, sharex=True)
ax2W.scatter(pTParams[1], Wmass_pred[1])
ax2W.errorbar(pTParams[1], Wmass_pred[1], yerr=Wmass_err[1], fmt='o')
ax2W.set_ylabel('Predicted mass of W (GeV)')

max2 = max(pTParams[1]) + 0.1*max(pTParams[1])
min2 = -0.1*max(pTParams[1])
ax2W.set_xlim([min2,max2])

ax2Chi.scatter(pTParams[1], chi2min_fit[1])
ax2Chi.set_ylabel('Test stat. of W template fit')
ax2Chi.set_xlabel('Mean percent pT smear (Gaussian)')

fig2.savefig("./plots/gaus_pTdependent_smear.png")
