import numpy as np
import ROOT as R
import csv
from decimal import *
import matplotlib.pyplot as plt

filename = "./rootfiles/build_compare_templatesWm.root"
input = R.TFile.Open(filename, "UPDATE")

nominalH = input.Get("Wm80.40")
targetH = input.Get("GausSmearWm0.018000")

nbinsN = nominalH.GetNbinsX()
nbinsT = targetH.GetNbinsX()
countsN = []
countsT = []
ratioTN = []
nomEntries = nominalH.GetEntries() - (nominalH.GetBinContent(0) + nominalH.GetBinContent(nbinsN)) 
tarEntries = targetH.GetEntries() - (targetH.GetBinContent(0) + targetH.GetBinContent(nbinsT)) 

print(nomEntries)
print(tarEntries)
print(targetH.GetEntries())
print(targetH.GetBinContent(0))
print(targetH.GetBinContent(nbinsT))

AdjustEntries = 0
ratiosum = 0
entrycountN = 0
entrycountT = 0

print(type(nominalH))
print(nominalH.GetBinContent(1))
print(nominalH.GetBinContent(3))
print(nominalH.GetBinContent(20))
print(nominalH.GetBinContent(40))
print(nominalH.GetNbinsX())

print("LOOPING OVER BINS")
for bin in range (1, nbinsN + 1):
    
    bincountN = nominalH.GetBinContent(bin)
    countsN.append(bincountN)
    entrycountN += bincountN 
    
    bincountT = targetH.GetBinContent(bin)
    #print('target bin ' + str(bin) + 'content: ' + str(bincountT)
    bincountnormed = bincountT*(nomEntries/tarEntries)
    #print('after normalization: ' +  str(bincountnormed))
    countsT.append(bincountnormed)
    entrycountT += bincountnormed

    countratio = bincountnormed/bincountN
    ratioTN.append(countratio)
    ratiosum += countratio - 1
    

#xlims = nominalH.GetXaxis().GetXbins()
#print(type(xlims))
print('RATIO SUM: ' + str(ratiosum))
print(entrycountN)
print(entrycountT)

# Append one extra item to match the x axis created by np.linspace
countsN.append(1)
countsT.append(1)
ratioTN.append(1)

histx = np.linspace(30,50,nbinsN+1, endpoint=True)
normline = [1 for i in range(nbinsN+1)]


print(len(histx))
print(histx[-1])
print(histx[0])
print(histx[1])
print(histx[2])
#plt.bar(histx, bincounts, 0.5, align='edge', fill=False)

#plt.show()

fig1, (axHist,axRatios) = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[5,3]})

axHist.step(histx, countsN, where='post')
axHist.step(histx, countsT, where='post')
axHist.set_xlim(30,50)
axHist.set_ylim(0,max(countsN)*1.1)
axHist.set_ylabel('Event Count')

axRatios.plot(histx,normline)
axRatios.step(histx, ratioTN, where='post')
axRatios.set_xlabel('mu_PT')
axRatios.set_ylabel('Target/Nominal')

fig1.savefig('plots/compare_hists.png')

#for bin in range (:
 #   bin = nominalH.GetBinContent(i)

#print(len(bincounts))
#print(bincounts[14]) 
    

