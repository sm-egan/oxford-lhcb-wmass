#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TCanvas.h"
#include "math.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include "TSystem.h"
#include "TH2F.h"
#include <vector>
#include "TemplateStruct.h"


using namespace std;

int main(){
  int hist_dims[3] = {40,30,50};
  int ndata_files = 3;
  int ntemplates = 11;
  Double_t Mnom = 80.4;
  Double_t gamma = 2.15553;
  
  vector<TH1F *> toys, templates;
  vector<Double_t> Wmasses, pTparams;

  vector<TH1F *>::iterator toyit, templateit;
  vector<Double_t>::iterator Wmassit, pTparamsit;

  string output_name = "~/oxford-lhcb-wmass/rootfiles/build_compare_templates.root";

  TH1::SetDefaultSumw2();

  //define histogram
  TH1F *h_muPT= new TH1F("h_mu_PT","#mu P_{T}",hist_dims[0], hist_dims[1], hist_dims[2]);

  //create output file  	
  TFile *output = new TFile(output_name.c_str(),"RECREATE");
 
  // load root files 

  TChain *MCDecayTree = new TChain("MCDecayTree");
  char filename[200];

  for(int k=1;k<=ndata_files;k++){
    if (k<10){
      sprintf(filename,"/data/lhcb/users/pili/forShannon/13000__Wp__PowhegPythia__as0.138_IKT1.0__evts0__seed000%u.root",k);
    } else sprintf(filename,"/data/lhcb/users/pili/forShannon/13000__Wp__PowhegPythia__as0.138_IKT1.0__evts0__seed00%u.root",k);
    cout <<"filename is  " << filename << endl;
    MCDecayTree->Add(filename);
  }
  
//Declaration of leaves types
  Float_t prop_M, mu_PT;
  MCDecayTree->SetBranchAddress("mu_PT", &mu_PT);
  MCDecayTree->SetBranchAddress("prop_M", &prop_M);


  stringstream snominal, sreweight;
  snominal << fixed << setprecision(2) << Mnom;
  string template_name;
  string toy_name = "pTsmear";

  cout << "error marker 1" << endl;

  // Initialize a vector of W mass hypotheses which will assist in iterating over and filling template histograms
  for (Double_t Mhyp=79.8; Mhyp<=80.8; Mhyp+=(80.8-79.8)/(ntemplates-1)) { //Mhyp for mass hypothesis
    Wmasses.push_back(Mhyp);  
    sreweight << fixed << setprecision(2) << Mhyp;
    template_name = "Reweight"+ sreweight.str() +"Nominal"+ snominal.str();

    TH1F *hweighted_template = new TH1F(template_name.c_str(), "mu_PT", hist_dims[0], hist_dims[1], hist_dims[2]);    
    templates.push_back(hweighted_template);
    
    //Resetting the stringstream
    sreweight.str(string());
  }

  for (Double_t pTparam = 1.0; pTparam <=1.0; pTparam += 0.5) {
    pTparams.push_back(pTparam);
    TH1F *hpT_toy = new TH1F(toy_name.c_str(), "mu_PT", hist_dims[0], hist_dims[1], hist_dims[2]);
    toys.push_back(hpT_toy);
    toy_name = "";
  }

  Long64_t nentries = MCDecayTree->GetEntries();
  cout << "nentries: " << nentries << endl;
  Long64_t nbytes = 0;
  
  int templateindex = 0;
  int toyindex = 0;

  //  cout << "Error marker 2" << endl;

  for (Long64_t eventit=0; eventit < 1000; eventit++) { //usually i< nentries*split_ratio for full data set
    nbytes += MCDecayTree->GetEntry(eventit);
    
    //cout << "error marker 3" << endl;
    cout << "currently at " << eventit << "th iteration" << endl;

    if(mu_PT > 30 && mu_PT < 50){ // apply a cut
      if (eventit%2 == 0){
	for (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	  //cout << "error marker 4" <<endl;
	  toys[0]->Fill(mu_PT);
	}
      } else {
	for (Wmassit = Wmasses.begin(); Wmassit != Wmasses.end(); ++Wmassit) {
	  //cout << "error marker 5" <<endl;
	  templates[templateindex]->Fill(mu_PT, TMath::BreitWigner(mu_PT, *Wmassit, gamma)/TMath::BreitWigner(mu_PT, Mnom, gamma));
	  templateindex = (templateindex + 1) % templates.size();
	  //cout << "templateindex is currently: " << templateindex <<endl;
	}
      }
    }
  }
 
  cout << "Error marker 6" << endl;
  
  //plot_overlaid_histograms();

  output->Write();
  output->Close();
  return 0;
}
