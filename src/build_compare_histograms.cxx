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

void plot_overlaid_histograms (vector<TH1F *> template_vect, TFile output) {
  bool firsthist = true;
  int colourit = 1;
  vector<TH1F *>::iterator templateit;
    // Below might be necessary for larger number of hypotheses
    //int lineit = 1;
    
  TCanvas* c = new TCanvas("cBW", "mu PT with different W mass hypotheses");
      
  for (templateit = this->template_vect.begin(); templateit != this->template_vect.end(); templateit++) {
      //brackets around *templateit ensure that we are acting on the pointer to the TH1F, not the iterative pointer to the pointer 
    (*templateit)->SetLineColor(colourit);

    if(firsthist) {
      (*templateit)->Draw();
      firsthist = false;
    } else (*templateit)->Draw("SAME");
  }
  colourit++;

  c->Print("~/oxford-lhcb-wmass/plots/templatesHist.png");
  c->Print("~/oxford-lhcb-wmass/plots/templatesHist.pdf");
  c->Close();    
}

void template_chi2 (vector<TH1F *> template_vect, vector<TH1F *> toy_vect, vector<Double_t> Wmass_vect, TFile output) {
  
  vector<TH1F *>::iterator templateit;
  vector<TH1F *>::iterator toyit;
  int ntemplates = template_vect.size();

  cout << "All vectors and iterators initialized.  template_vect is of length " << template_vect.size() << ". toy_vect is of length " <<  toy_vect.size() << ". Wmass_vect is of length " << (this->Wmasses).size() << endl;
  cout << template_vect.front() << endl;
  char scaledhname[200];

  int toyindex=0, templateindex=0;  
  // Copy template hypothesis info onto an array of fixed size to hopefully be accepted by TGraph->DrawGraph()
  
  Double_t Wmass_arr[ntemplates];
  cout << "Copying TemplateStruct->Wmasses to an array" << endl;
  
  TGraph *chi2Plot = new TGraph(ntemplates);
  //char chi2plot_name[1000];
  string chi2plot_name = "chi2 plot";
  Double_t chi2point;
  /*
  Double_t par0, par1, par2;
  TF1 *chi2fit;
  char chi2fit_name[100];  
  */
//Unit normalize all of the toy vectors and perform chi square test, saving the data to an array as you go
  TCanvas* c2 = new TCanvas("cchi2", "chi2 with different W mass hypotheses");

  for (toyit = toy_vect.begin(); toyit != toy_vect.end(); toyit++) {
    cout << "performing loop over toys" << endl;
    if ((*toyit)->GetSumw2N() == 0) {
      cout << "Warning! Weights do not seem to be stored" << endl;
    }

    for (templateit = template_vect.begin(); templateit != template_vect.end(); templateit++){

      (*toyit)->Scale((*templateit)->Integral()/(*toyit)->Integral());
      sprintf(scaledhname, "%s_scaledto_%s", (*toyit)->GetName(), (*templateit)->GetName());
      output->WriteTObject(*toyit, scaledhname, "Overwrite");
      chi2point = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
      
      cout << "Copying chi2 result for toy " << toyindex << " and template " << templateindex << ": " << chi2point << endl;
      chi2Plot->SetPoint(templateindex,Wmass_vect[templateindex],chi2point);

      templateindex = (templateindex + 1) % template_vect.size();
    }
    cout << "Drawing graph of toy " << toyindex << endl;
    //chi2Plot->DrawGraph(template_vect.size(), Wmass_arr, chi2_results[toyindex]);
    chi2Plot->SetMarkerColor(4);
    chi2Plot->SetMarkerStyle(8);
    chi2Plot->Draw("AP");
    
    
    cout << "Writing chi2 plot to root file." << endl;
    //sprintf(chi2plot_name, "chi2_plot_%s_vs_%s", (*toyit)->GetName(), (*templateit)->GetName());
    //cout << "Name of chi2 plot has been assigned" << endl;
    output->WriteTObject(chi2Plot, chi2plot_name.c_str(), "Overwrite");
    cout << "Write successful" << endl;
    /*
    sprintf(chi2fit_name, "chi2fit_%s_vs_%s", (*toyif)->GetName(), (*templateit)->GetName());
    chi2fit = new TF1(chi2fit_name, "[0]*x*x+[1]*x+[2]");
    chi2Plot->Fit(chi2fint_name, "", "", 79.9, 80.8);
    */
    toyindex++;
  }
  cout << "Saving chi2 plot as an image." << endl;
  c2->Print("./plots/chi2Plot.png");
  c2->Close();  
}


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
  }

  for (Double_t pTparam = 1.0; pTparam <=1.0; pTparam += 0.5) {
    pTparams.push_back(pTparam);
    TH1F *hpT_toy = new TH1F(toy_name.c_str(), "mu_PT", hist_dims[0], hist_dims[1], hist_dims[2]);
    toys.push_back(hpT_toy);
  }

  Long64_t nentries = MCDecayTree->GetEntries();
  Long64_t nbytes = 0;
  
  int templateindex = 0;
  int toyindex = 0;

  cout << "Error marker 2" << endl;

  for (Long64_t eventit=0; eventit < nentries; eventit++) { //usually i< nentries*split_ratio for full data set
    nbytes += MCDecayTree->GetEntry(i);
    
    cout << "error marker 3" << endl;
    cout << "currently at " << i << "th iteration" << endl;

    if(mu_PT > 30 && mu_PT < 50){ // apply a cut
      if (eventit%2 == 0){
	for (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	  cout << "error marker 4" <<endl;
	  toys[0]->Fill(mu_PT);
	}
      } else {
	for (Wmassit = Wmasses.begin(); Wmassit != Wmasses.end(); ++Wmassit) {
	  cout << "error marker 5" <<endl;
	  templates[templateindex]->Fill(mu_PT, TMath::BreitWigner(mu_PT, *Wmassit, gamma)/TMath::BreitWigner(mu_PT, Mnom, gamma));
	  templateindex = (templateindex + 1) % templates.size();
	  cout << "templateindex is currently: " << templateindex <<endl;
	}
      }
    }
  }
 
  cout << "Error marker 4" << endl;
  
  plot_overlaid_histograms();

  output->Write();
  output->Close();
  return 0;
}
