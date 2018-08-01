#from ROOT import gROOT as R
#R.SetBatch(True) 
import sys
import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt
from plot_testing import hist_ratio_plot

#### MANUALLY SET THE DESIRED FILE AND NUMBER OF PT ADJUSTMENT METHODS HERE ###
fileinfo = sys.argv[1]
filename = './rootfiles/' + fileinfo + '.root'
#pTmethods = ['GausSmear']
pTmethods = ['GausSmear', 'GausSmear_pTdependent', 'ConstFactor', 'CurveOffset']
#pTmethods = ['GausSmear', 'GausSmear_pTdependent', 'ConstFactor']

if filename.find("Wm") > -1:
    Wcharge = "Wm"
elif filename.find("Wp") > -1:
    Wcharge = "Wp"
else:
    print("Could not find a keyword for W charge in filename")

npTmethods = len(pTmethods)


input = R.TFile.Open(filename, "UPDATE")
pTParams = [[]]
Wmass_pred = [[]*npTmethods]
Wmass_err = [[]*npTmethods]
chi2min_fit = [[]*npTmethods]

nominalH_name = Wcharge + "80.40"

### INITIALIZE THE FIT FUNCTION AND NAME THE PARAMETERS ACCORDINGLY ###
quadfit = R.TF1("chi2_quadfit", "[0]+(1/([2]**2))*(x-[1])**2", 79.8, 80.8)
quadfit.SetParName(0,"MinChi2")
quadfit.SetParName(1,"WMass")
quadfit.SetParName(2,"WMassUncertainty")


### READ IN THE PT PARAMETERS FROM CSV FILE OUTPUT BY build_compare_histograms.exe ###
firstrow = True
parameterfile = './pTparameters/' + fileinfo + '.csv'  
with open(parameterfile) as csvfile:
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
        
        chi2Plot_name = "chi2plot" + Wcharge + str(toyset) + str(toyit)    
        targetH_name = pTmethods[toyset] + Wcharge + str(toyit)
        print("Loading graphs: " + chi2Plot_name + " and " + targetH_name)
        
        print(type(input.Get(targetH_name)))
        
        chi2Plot = input.Get(chi2Plot_name)
        chi2stats = chi2Plot.GetY()
    
        p0_init = min(chi2stats)
        print("p0_init: " + str(p0_init))
        p1_init = chi2Plot.GetX()[np.argmin(chi2stats)]
        print("p1_init: " + str(p1_init))
        p2_init = 0.01

        wH_name = Wcharge + "%.2f" % p1_init
        print("ADDING HIST COMPARISON OF " + wH_name)
        nominalH_name = Wcharge + "80.40"
        hist_ratio_plot(input, nominalH_name, targetH_name, wH_name)    

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

        chi2imagepdf = "./plots/chi2fits/chi2Plot" + Wcharge + str(toyset) +  str(toyit) + "_fit.pdf"
        chi2imagepng = "./plots/chi2fits/chi2Plot" + Wcharge + str(toyset) +  str(toyit) + "_fit.png"
        c.Print(chi2imagepdf)
        c.Print(chi2imagepng)
        c.Close()

    '''
    if not ((len(Wmass_pred) == len(Wmass_err)) and (len(Wmass_pred) == len(chi2min_fit)) and (len(Wmass_pred) == len(pTParams))):
        print("Warning!  The arrays you wish to plot are not of equal length!")
        print('Wmass predictions: ' + str(len(Wmass_pred)))
        print('Wmass errors: ' + str(len(Wmass_err)))
        print('pT parameters: ' + str(len(pTParams)))
        print('chi2 test statistics: ' + str(len(chi2min_fit)))
    '''
    print('Plotting W predictions and minimum chi square for pT method ' + str(toyset))
    fig, (axW, axChi) = plt.subplots(2, sharex=True)
    axW.scatter(pTParams[toyset], Wmass_pred[toyset])
    axW.errorbar(pTParams[toyset], Wmass_pred[toyset], yerr=Wmass_err[toyset], fmt='o')
    axW.set_ylabel('Predicted mass of W (GeV)')
    
    setmax = max(pTParams[toyset])
    setmin = min(pTParams[toyset])

    xmax = setmax + (setmax-setmin)*0.1
    xmin = setmin-(setmax-setmin)*0.1
    axChi.set_xlim([xmin,xmax])

    axChi.scatter(pTParams[toyset], chi2min_fit[toyset])
    
    axChi.set_ylabel(r'$\chi ^2$ Test stat. of fit')
    xlabelstr = pTmethods[toyset] + ' parameter'
    axChi.set_xlabel(xlabelstr)

    fig.suptitle(Wcharge + ' samples')
    if toyset == 0:
        axW.set_title(r"$\frac{\sigma p}{p} \propto x$")
    if toyset == 1:
        axW.set_title(r"$\frac{\sigma p}{p} \propto x*p$")
    if toyset == 2:
        axW.set_title(r"$p \rightarrow x*p$")
    if toyset == 3:
        axW.set_title(r"$\frac{q}{p} \rightarrow \frac{q}{p} + x$")
    fig.savefig('./plots/Wpredictionplots/' + Wcharge + 'predictions_pTmethod' + str(toyset) + '.pdf', bbox_inches='tight')
    fig.savefig('./plots/Wpredictionplots/' + Wcharge + 'predictions_pTmethod' + str(toyset) + '.png', bbox_inches='tight')

    plt.close('all')
