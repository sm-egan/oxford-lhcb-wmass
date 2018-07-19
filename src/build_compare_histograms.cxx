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
  TemplateStruct *ts = new TemplateStruct();

  TH1::SetDefaultSumw2();

  //define histogram
  TH1F *h_muPT= new TH1F("h_mu_PT","#mu P_{T}",hist_dims[0], hist_dims[1], hist_dims[2]);

  //create output file  	
  TFile *output = new TFile(ts->output_name,"RECREATE");
 
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

  long int maxentries = lrint((ts->split_ratio)*(MCDecayTree->GetEntries()));

  //for (loop over particular toy aspects) { 
  MCDecayTree->Draw("mu_PT >> h_mu_PT", "(mu_PT > 30) && (mu_PT < 50)", "", maxentries);
  ts->toys.push_back(h_muPT);
  //} 
  
  for (Double_t Mhyp=79.8; Mhyp<=80.8; Mhyp+=(80.8-79.8)/(ntemplates-1)) { //Mhyp for mass hypothesis
      ts->drawH_BW(MCDecayTree, "prop_M", "mu_PT", hist_dims, Mnom, Mhyp);
      
      /* The goal of the lines commented out here was to make sure that there would always be n templates, even when the increment of the addition does not divide perfectly into 1
      if (((Mhyp +(80.8-79.8)/(ntemplates-1)) > 80.8) && (Mhyp != 80.8)) {
	this->drawH_BW(MCDecayTree, "prop_M", "mu_PT", hist_dims, Mnom, 80.8);
	}*/
  }
  output->Write();
  ts->template_chi2();
  
  //  ts->create_templates(hist_dims, 11, 2);
  cout << "template_vect is of length " << ts->templates.size() << ". toy_vect is of length " <<  ts->toys.size() << ". Wmass_vect is of length " << ts->Wmasses.size() << endl;
  cout << (ts->templates)[0] << endl;
  cout << ts->toys.front() << endl;
  //cout << "templates created. Now performing chi2 test" << endl;
  //ts->template_chi2();
}
