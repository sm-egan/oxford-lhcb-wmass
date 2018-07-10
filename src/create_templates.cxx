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




using namespace std;

TH1F* BWweight (TChain* EventChain, string ReweightBranch, string HistBranch, int hist_dims[3], Double_t nominal_mean, Double_t reweight_mean, Double_t gamma) {
  //Problem, the EventTree object passed is a TTree and not a TChain, this function may therefore inflexible for adding information over many trees
  string s = "Reweight"+to_string(reweight_mean)+"Nominal"+to_string(nominal_mean);
  //int n = s.length();
  //char char_array[n+1];
  //strcpy(char_array, s.c_str());

  TH1F *hweighted = new TH1F(s.c_str(), HistBranch.c_str(), hist_dims[0], hist_dims[1], hist_dims[2]);
  EventChain->Draw("EventChain->GetBranch(HistBranch)>>hweighted",
		   "(HistBranch > hist_dims[1] && HistBranch < hist_dims[2])*(TMath::BreitWigner(ReweightBranch, reweight_mean, gamma)/TMath::BreitWigner(TBranch, nominal_mean, gamma))"); //It's probably the math expressions here with HistBranch which cause it to fail
  return hweighted;
}

void create_templates(){

  TH1::SetDefaultSumw2();

  //define histograms
  int hist_dims[3] = {40,30,50};
  TH1F *h_muPT= new TH1F("h_mu_PT","#mu P_{T}",hist_dims[0], hist_dims[1], hist_dims[2]);
 
  //create output file
  TString name= "/home/egan/oxford-lhcb-wmass/rootfiles/create_templates.root";  	
  TFile *output = new TFile(name,"RECREATE");
 
  // load root files 

  TChain *MCDecayTree = new TChain("MCDecayTree");
  char filename[200];

  for(int k=1;k<3;k++){
    sprintf(filename,"/data/lhcb/users/pili/forShannon/13000__Wp__PowhegPythia__as0.138_IKT1.0__evts0__seed000%i.root",k);
      cout <<"filename is  " << filename << endl;
      MCDecayTree->Add(filename);}
  
      /* for(int k=10;k<21;k++){
    sprintf(filename,"/data/lhcb/users/pili/forShannon/13000__Wp__PowhegPythia__as0.138_IKT1.0__evts0__seed00%i.root",k);
    cout <<"filename is " << filename << endl;
    MCDecayTree->Add(filename);
    }*/

//Declaration of leaves types
   Float_t         prop_M;
   Float_t         dilep_PT;
   Float_t         dilep_M;
   Float_t         dilep_Y;
   Float_t         mu_PT;
   Float_t         mu_ETA;
   Float_t         mu_PHI;
   Int_t           mu_ID;
   Float_t         nu_PT;
   Float_t         nu_ETA;
   Float_t         nu_PHI;
   Int_t           nu_ID;
   Int_t           id1pdf;
   Int_t           id2pdf;
   Float_t         x1pdf;
   Float_t         x2pdf;

   // Set branch addresses.
   MCDecayTree->SetBranchAddress("prop_M",&prop_M);
   MCDecayTree->SetBranchAddress("dilep_PT",&dilep_PT);
   MCDecayTree->SetBranchAddress("dilep_M",&dilep_M);
   MCDecayTree->SetBranchAddress("dilep_Y",&dilep_Y);
   MCDecayTree->SetBranchAddress("mu_PT",&mu_PT);
   MCDecayTree->SetBranchAddress("mu_ETA",&mu_ETA);
   MCDecayTree->SetBranchAddress("mu_PHI",&mu_PHI);
   MCDecayTree->SetBranchAddress("mu_ID",&mu_ID);
   MCDecayTree->SetBranchAddress("nu_PT",&nu_PT);
   MCDecayTree->SetBranchAddress("nu_ETA",&nu_ETA);
   MCDecayTree->SetBranchAddress("nu_PHI",&nu_PHI);
   MCDecayTree->SetBranchAddress("nu_ID",&nu_ID);
   MCDecayTree->SetBranchAddress("id1pdf",&id1pdf);
   MCDecayTree->SetBranchAddress("id2pdf",&id2pdf);
   MCDecayTree->SetBranchAddress("x1pdf",&x1pdf);
   MCDecayTree->SetBranchAddress("x2pdf",&x2pdf);

   Long64_t nentries = MCDecayTree->GetEntries();
   
   Double_t gamma = 2.15553;
   Double_t Mnom = 80.3819;
   
   Long64_t nbytes = 0;
   for (Long64_t i=0; i<(nentries/10);i++) { //i usually nentries for full data set

    nbytes += MCDecayTree->GetEntry(i);
    
    if(mu_PT > 30 && mu_PT < 50){ // apply a cut
      h_muPT->Fill(mu_PT);} // add data from each n-tuple to the same histogram
	// h_muPT->Fill(mu_PT,weight); if you want to give a weight to the histogram 
  }


   TCanvas *c0 = new TCanvas("c0","c0");
   c0->cd();
   
     h_muPT->Draw();
     output->WriteTObject(h_muPT,h_muPT->GetName(),"Overwrite");
     BWweight(MCDecayTree, "prop_M", "mu_PT", hist_dims, 80.4, 80.6, gamma);

  c0->Close();   
  output->Close();
   
}

int main(){
  create_templates();
}
