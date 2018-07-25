#from ROOT import gROOT as R
#R.SetBatch(True) 
import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt


#### MANUALLY SET THE DESIRED FILE AND NUMBER OF PT ADJUSTMENT METHODS HERE ###
filename = "./rootfiles/build_compare_templatesWm.root"
input = R.TFile.Open(filename, "UPDATE")
npTmethods = 2


pTParams = [[]]
Wmass_pred = [[]*npTmethods]
Wmass_err = [[]*npTmethods]
chi2min_fit = [[]*npTmethods]

Wcharge = "Wm"

### INITIALIZE THE FIT FUNCTION AND NAME THE PARAMETERS ACCORDINGLY ###
quadfit = R.TF1("chi2_quadfit", "[0]+(1/([2]**2))*(x-[1])**2", 79.8, 80.8)
quadfit.SetParName(0,"MinChi2")
quadfit.SetParName(1,"WMass")
quadfit.SetParName(2,"WMassUncertainty")


### READ IN THE PT PARAMETERS FROM CSV FILE OUTPUT BY build_compare_histograms.exe ###
firstrow = True 
with open('./pTparameters.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC)
    for row in readCSV:
        if (firstrow): 
            pTParams[0] = row
            firstrow = False
        else:
            pTParams.append(row)

print('Imported pT parameters')
print(pTParams)


### LOOP OVER THE PT ADJUSTMENT METHODS AND PRODUCE A 2-SUBPLOT FIGURE SHOWING BEST FIT W MASS PREDICTION AND MINIMUM CHI SQUARE FOR EACH ###
for toyset in range (0,npTmethods):
    print('Now plotting set of toys: ' + str(toyset))
    ntoys = len(pTParams[toyset])
    for toyit in range (0, ntoys):    
        print('Now plotting toy' + str(toyit))
        if filename.find("Wm") > -1:
            chi2Plot_name = "chi2plot" + Wcharge + str(toyset) + str(toyit)    
            print("Loading graph: " + chi2Plot_name)
        else:
            Wcharge = "Wp"
            chi2Plot_name = "chi2plot" + Wcharge + str(toyset) + str(toyit)    
            
        chi2Plot = input.Get(chi2Plot_name)
        chi2stats = chi2Plot.GetY()
    
        p0_init = min(chi2stats)
        print("p0_init: " + str(p0_init))
        p1_init = chi2Plot.GetX()[np.argmin(chi2stats)]
        print("p1_init: " + str(p1_init))
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

        chi2image = "./plots/chi2Plot" + Wcharge + str(toyset) +  str(toyit) + "_fit.png"
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
    print('Plotting W predictions and minimum chi square for toy ' + str(toyset))
    fig, (axW, axChi) = plt.subplots(2, sharex=True)
    axW.scatter(pTParams[toyset], Wmass_pred[toyset])
    axW.errorbar(pTParams[toyset], Wmass_pred[toyset], yerr=Wmass_err[toyset], fmt='o')
    axW.set_ylabel('Predicted mass of W (GeV)')
    
    setmax = max(pTParams[toyset])

    xmax = setmax + 0.1*setmax
    xmin = -0.1*setmax
    axW.set_xlim([xmin,xmax])

    axChi.scatter(pTParams[toyset], chi2min_fit[toyset])
    axChi.set_ylabel('Test stat. of W template fit')
    axChi.set_xlabel('pT adjustment parameter')

    fig.savefig('./plots/' + Wcharge + 'predictions_pTmethod' + str(toyset) + '.png')



### PLOT HISTOGRAMS OF TOYS OVERLAID WITH NOMINAL WMASS HISTOGRAM, ALONG WITH RATIOS ###
