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


TemplateStruct::TemplateStruct () {
	this->split_ratio = 2;	
	this->output_name = "~/oxford-lhcb-wmass/rootfiles/create_templates.root";
}

TemplateStruct::TemplateStruct (Double_t split, TString root_file) {
	this->split_ratio = split;
	this->output_name = root_file;
}

void TemplateStruct::drawH_BW (TChain* EventChain, string ReweightBranch, string HistBranch, int hist_dims[3], Double_t nominal_mean, Double_t reweight_mean) {
  //Set the title of the histogram to include the Wmass hypothesis
  stringstream snominal, sreweight;
  snominal << fixed << setprecision(2) << nominal_mean;
  sreweight << fixed << setprecision(2) << reweight_mean;

  // TRUNCATE DECIMAL REALLY CAN'T BE METHOD OF TEMPLATESTRUCT, PROBABLY HAVE TO DO THIS MANUALLY
  string s = "Reweight"+ sreweight.str() +"Nominal"+ snominal.str();
  Double_t gamma = 2.15553;
  //Long64_t nentries = EventChain->GetEntries();

  TH1F *hweighted = new TH1F(s.c_str(), HistBranch.c_str(), hist_dims[0], hist_dims[1], hist_dims[2]);
  
  //Create strings of characters to send to the draw expression which integrate the passed parameters
  char varexp[100];
  sprintf(varexp, "%s >> %s", HistBranch.c_str(), s.c_str());
  char selection[300];
  sprintf(selection, "%s > %u && %s < %u", HistBranch.c_str(), hist_dims[1], HistBranch.c_str(), hist_dims[2]);
  char weightexp[300];
  sprintf(weightexp, "(%s)*(TMath::BreitWigner(%s, %f, %f)/TMath::BreitWigner(%s,%f,%f))", 
	  selection,ReweightBranch.c_str(),reweight_mean,gamma,ReweightBranch.c_str(),nominal_mean,gamma);

  EventChain->Draw(varexp,weightexp); 
  this->templates.push_back(hweighted);
  this->Wmasses.push_back(reweight_mean);
}

void TemplateStruct::create_templates(int hist_dims[3], int ntemplates=5, int ndata_files=2) {

  TH1::SetDefaultSumw2();

  //define histogram
  TH1F *h_muPT= new TH1F("h_mu_PT","#mu P_{T}",hist_dims[0], hist_dims[1], hist_dims[2]);

  //create output file
  TString name= "~/oxford-lhcb-wmass/rootfiles/create_templates.root";  	
  TFile *output = new TFile(name,"RECREATE");
 
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
   Double_t Mnom = 80.3819;

   Long64_t nbytes = 0;
   //vector<TH1F *> hist_vect; Not needed in class implementation since it is given from the beginning

   for (Long64_t i=0; i<(nentries/split_ratio);i++) { //usually i< nentries for full data set

    nbytes += MCDecayTree->GetEntry(i);
    
    if(mu_PT > 30 && mu_PT < 50){ // apply a cut
      h_muPT->Fill(mu_PT);} // add data from each n-tuple to the same histogram
	// h_muPT->Fill(mu_PT,weight); if you want to give a weight to the histogram 
  }
   
    h_muPT->Draw();
    output->WriteTObject(h_muPT,h_muPT->GetName(),"Overwrite");
    //vector<Double_t> Wmass_vect;

    for (Double_t Mhyp=79.8; Mhyp<=80.8; Mhyp+=(80.8-79.8)/(ntemplates-1)) { //Mhyp for mass hypothesis
      this->drawH_BW(MCDecayTree, "prop_M", "mu_PT", hist_dims, Mnom, Mhyp);
    }

    bool firsthist = true;
    int colourit = 1;
    // Below might be necessary for larger number of hypotheses
    //int lineit = 1;
    
    TCanvas* c = new TCanvas("cBW", "mu PT with different W mass hypotheses");
      
    for (vector<TH1F *>::iterator histit = this->templates.begin(); histit != this->templates.end(); histit++, colourit++) {
      //brackets around *histit ensure that we are acting on the pointer to the TH1F, not the iterative pointer to the pointer  
      output->WriteTObject(*histit, (*histit)->GetName(),"Overwrite");
      (*histit)->SetLineColor(colourit);

      if(firsthist) {
	(*histit)->Draw();
	firsthist = false;
      } else (*histit)->Draw("SAME");
    }
    c->Print("~/oxford-lhcb-wmass/plots/WmasshypHist.png");
    c->Print("~/oxford-lhcb-wmass/plots/WmasshypHist.pdf");
    c->Close();
    
    //TemplateStruct ts; This function is already acting on an instance
    //ts.Wmasses = Wmass_vect; Should be filled already in the drawH_BW function
    //ts.templates = hist_vect; "      
    //ts.toys.push_back(h_muPT); 
    this->toys.push_back(h_muPT);
}