#from ROOT import gROOT as R
#R.SetBatch(True) 
import ROOT as R
input = R.TFile.Open("./rootfiles/create_templates.root", "UPDATE")

chi2Plot = input.Get("chi2 plot")

quadfit = R.TF1("chi2_quadfit", "[0]+[2]*(x-[1])**2", 79.8, 80.8)
quadfit.SetParameters(33,80.4,10000)

chi2Plot.Fit("chi2_quadfit")

chi2min = quadfit.GetMinimum()
Mw = quadfit.GetX(chi2min)
MwSigma = abs(quadfit.GetX(chi2min+1)-Mw)

print('The best W mass is ' + str(Mw) + '+/-' + str(MwSigma))

c = R.TCanvas()
chi2Plot.Draw("AP")

c.Print("./plots/chi2Plot.png")
c.Close()
