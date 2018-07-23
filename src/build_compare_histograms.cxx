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
#include "TRandom.h"
#include <fstream>

using namespace std;

/*
void write_vector (ofstream bufferstream, vector<Double_t> data_vect) {

  if (!bufferstream.is_open()) {
      cout << "ERROR: ofstream not associated to open file" << endl;  
  } else {
    for (vector<Double_t>::iterator data_it = data_vect.begin(); data_it != data_vect.end(); ++data_it) { 
      if (data_it == data_vect.begin()) {
	bufferstream << *data_it; 
      } else {
	bufferstream << ',' << *data_it; 
      }
    }
    cout << '\n';
  }
}
*/ 

int main(){
  int hist_dims[3] = {40,30,50};
  int ndata_files = 1;
  int ntemplates = 11;
  int ntoys = 6;
  Double_t Mnom = 80.4;
  Double_t gamma = 2.15553;
  
  vector<TH1F *> templates;
  vector< vector<TH1F *> > toys(2);
  vector<Double_t> Wmasses;
  vector< vector<Double_t> > pTparams(3);

  vector<TH1F *>::iterator toyit, templateit;
  vector<Double_t>::iterator Wmassit, pTparamsit;

  string output_name = "~/oxford-lhcb-wmass/rootfiles/build_compare_templates.root";

  TH1::SetDefaultSumw2();

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

  MCDecayTree->SetBranchStatus("*", 0);
  MCDecayTree->SetBranchStatus("mu_PT", 1);
  MCDecayTree->SetBranchStatus("prop_M",1);

  stringstream snominal, sreweight;
  snominal << fixed << setprecision(2) << Mnom;
  string template_name;
  string toy_name1 = "GausSmear";
  string toy_name2 = "GausSmear_pTdependent";

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

  cout << "error marker 2" << endl;
  
  for (Double_t pTparam = 0.0; pTparam <=0.03; pTparam += (0.03)/(ntoys-1)) {
    pTparams[0].push_back(pTparam);
    toy_name1 = toy_name1 + to_string(pTparam);
    TH1F *hpT_toy = new TH1F(toy_name1.c_str(), "mu_PT with simple momentum smear", hist_dims[0], hist_dims[1], hist_dims[2]);
    toys[0].push_back(hpT_toy);
    toy_name1 = "GausSmear";
  }

  cout << "error marker 3" << endl;

  for (Double_t pTparam = 0.0; pTparam <= 0.00075; pTparam += (0.00075)/(ntoys-1)) {
      pTparams[1].push_back(pTparam);
      toy_name2 = toy_name2 + to_string(pTparam);
      TH1F *hpT_toy = new TH1F(toy_name2.c_str(), "mu_PT with momentum-dependent momentum smear", hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[1].push_back(hpT_toy);
      toy_name2 = "GausSmear_pTdependent";
  }

  
  ofstream pT_ofs ("./pTparameters.csv", ofstream::out);
  cout << "error marker 4" << endl;

  for (int vect_row = 0; vect_row < pTparams.size(); ++vect_row) {
    for (vector<Double_t>::iterator data_it = pTparams[vect_row].begin(); data_it != pTparams[vect_row].end(); ++data_it) { 
      if (data_it == pTparams[vect_row].begin()) {
	pT_ofs << *data_it; 
      } else {
	pT_ofs << ',' << *data_it; 
      }
    }
    pT_ofs << endl;  
  }
  
  pT_ofs.close();

  Long64_t nentries = MCDecayTree->GetEntries();
  cout << "nentries: " << nentries << endl;
  Long64_t nbytes = 0;
  int toyindex =0, templateindex=0;

  //  cout << "Error marker 2" << endl;

  for (Long64_t eventit=0; eventit < nentries; eventit++) { //usually i< nentries*split_ratio for full data set
    nbytes += MCDecayTree->GetEntry(eventit);
    
    //cout << "error marker 3" << endl;
    //cout << "currently at " << eventit << "th iteration" << endl;

    //if(mu_PT > 30 && mu_PT < 50){ // apply a cut
      if (eventit%2 == 0){
	toyindex=0;
	for ( auto pTparam0 : pTparams[0] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	//cout << "error marker 4" <<endl;
	  mu_PT = mu_PT*gRandom->Gaus(1, pTparam0);
	  toys[0][toyindex]->Fill(mu_PT);
	  ++toyindex;
	}
	toyindex=0;
	for ( auto pTparam1 : pTparams[1] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	//cout << "error marker 4" <<endl;
	  mu_PT = mu_PT*gRandom->Gaus(1, pTparam1*mu_PT);
	  toys[1][toyindex]->Fill(mu_PT);
	  ++toyindex;
	}
      } else {
	templateindex = 0;
	const auto denominator = TMath::BreitWigner(prop_M, Mnom, gamma); 
	for ( auto Wmass : Wmasses ){ //Wmassit = Wmasses.begin(); Wmassit != Wmasses.end(); ++Wmassit) {
	  //cout << "error marker 5" <<endl;
	  templates[templateindex]->Fill(mu_PT, TMath::BreitWigner(prop_M, Wmass, gamma)/denominator);
	  //templateindex = (templateindex + 1) % templates.size();
	  ++templateindex;
	  //cout << "templateindex is currently: " << templateindex <<endl;
	}
      }
      //}
  }
 

  bool firsthist = true;
  int colourit = 1;
    // Below might be necessary for larger number of hypotheses
    //int lineit = 1;
    
  TCanvas* c = new TCanvas("cBW", "mu PT with different W mass hypotheses");
      
  for (templateit = templates.begin(); templateit != templates.end(); templateit++, colourit++) {
    (*templateit)->SetLineColor(colourit);

    if(firsthist) {
      (*templateit)->Draw();
      firsthist = false;
    } else (*templateit)->Draw("SAME");
  }
  c->Print("~/oxford-lhcb-wmass/plots/templatesHist.png");
  c->Print("~/oxford-lhcb-wmass/plots/templatesHist.pdf");
  c->Close();    

       /*
  TCanvas* c1 = new TCanvas("cPTsmear", "mu PT with Gaussian smear");
      
  for (toyit = toys[0].begin(); toyit != toys[0].end(); toyit++, colourit++) {
    (*toyit)->SetLineColor(colourit);

    if(firsthist) {
      (*toyit)->Draw();
      firsthist = false;
    } else (*toyit)->Draw("SAME");
  }
  c1->Print("~/oxford-lhcb-wmass/plots/toys1Hist.png");
  c1->Print("~/oxford-lhcb-wmass/plots/toys1Hist.pdf");
  c1->Close();    
       */

  // PERFORMING A CHI SQUARE TEST

  cout << "All vectors and iterators initialized.  template_vect is of length " << templates.size() << ". toy_vect is of length " <<  toys.size() << ". Wmass_vect is of length " << Wmasses.size() << endl;
  cout << templates.front() << endl;
  //char scaledhname[300];

    // Copy template hypothesis info onto an array of fixed size to hopefully be accepted by TGraph->DrawGraph()
  
  //Double_t Wmass_arr[ntemplates];
  //cout << "Copying TemplateStruct->Wmasses to an array" << endl;
  //copy((Wmasses).begin(), (Wmasses).end(), Wmass_arr);
  
  TGraph *chi2Plot = new TGraph(ntemplates);
  
  char chi2plot_name[50];
  char legend_header[300];
  Double_t chi2point;
  toyindex=0;
  templateindex=0;  

  TCanvas* c2 = new TCanvas("cchi2", "chi2 with different W mass hypotheses");

  for (toyit = toys[0].begin(); toyit != toys[0].end(); toyit++) {
    
    cout << "performing loop over toys" << endl;
    if ((*toyit)->GetSumw2N() == 0) {
      cout << "Warning! Weights do not seem to be stored" << endl;
    }

    for (templateit = templates.begin(); templateit != templates.end(); templateit++){
      
      (*toyit)->Scale((*templateit)->Integral()/(*toyit)->Integral());
      sprintf(chi2plot_name, "chi2plot0%u", toyindex);
      
      //output->WriteTObject(*toyit, scaledhname, "Overwrite");
      chi2point = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
      
      cout << "Copying chi2 result for toy " << toyindex << " and template " << templateindex << ": " << chi2point << endl;
      chi2Plot->SetPoint(templateindex,Wmasses[templateindex],chi2point);

      templateindex = (templateindex + 1) % templates.size();
    }
    
    
    cout << "Drawing graph of toy " << toyindex << endl;
    chi2Plot->SetMarkerColor(4);
    chi2Plot->SetMarkerStyle(8);
    chi2Plot->Draw("AP");
    
    auto *legend = new TLegend();
    sprintf(legend_header, "Chi2: Gaussian pT smearing %f against W mass hypotheses", pTparams[0][toyindex]);
    legend->SetHeader(legend_header, "C");
    legend->Draw();
    
    cout << "Writing chi2 plot to root file." << endl;
    output->WriteTObject(chi2Plot, chi2plot_name, "Overwrite");
    cout << "Write successful" << endl;
    
    cout << "Saving chi2 plot as an image." << endl;
    c2->Print("./plots/chi2Plot.png");
    c2->Clear();  

    toyindex++;
  }
  c2->Close();

  toyindex=0, templateindex=0;

  TCanvas* c3 = new TCanvas("chi2", "toy: pT-dependent smearing, templates: W mass hypotheses");

  for (toyit = toys[1].begin(); toyit != toys[1].end(); toyit++) {
    
    cout << "performing loop over toys" << endl;
    if ((*toyit)->GetSumw2N() == 0) {
      cout << "Warning! Weights do not seem to be stored" << endl;
    }

    for (templateit = templates.begin(); templateit != templates.end(); templateit++){
      
      (*toyit)->Scale((*templateit)->Integral()/(*toyit)->Integral());
      sprintf(chi2plot_name, "chi2plot1%u", toyindex);
      
      //output->WriteTObject(*toyit, scaledhname, "Overwrite");
      chi2point = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
      
      cout << "Copying chi2 result for toy " << toyindex << " and template " << templateindex << ": " << chi2point << endl;
      chi2Plot->SetPoint(templateindex,Wmasses[templateindex],chi2point);

      templateindex = (templateindex + 1) % templates.size();
    }
    
    
    cout << "Drawing graph of toy " << toyindex << endl;
    chi2Plot->SetMarkerColor(4);
    chi2Plot->SetMarkerStyle(8);
    chi2Plot->Draw("AP");
    
    auto *legend = new TLegend();
    sprintf(legend_header, "Chi2: Gaussian pT smearing %f against W mass hypotheses", pTparams[0][toyindex]);
    legend->SetHeader(legend_header, "C");
    legend->Draw();
    
    cout << "Writing chi2 plot to root file." << endl;
    output->WriteTObject(chi2Plot, chi2plot_name, "Overwrite");
    cout << "Write successful" << endl;
    
    cout << "Saving chi2 plot as an image." << endl;
    c3->Print("./plots/chi2Plot.png");
    c3->Clear();  

    toyindex++;
  }
  c3->Close();

  
  output->Write();
  output->Close();
  return 0;
}


