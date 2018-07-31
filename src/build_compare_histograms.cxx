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
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include <TLorentzVector.h>


///data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root

Double_t GausSmear (Double_t value, Double_t adjparameter) {
  value = value*gRandom->Gaus(1, adjparameter);
  return value;
}

Double_t GausSmear_pTdependent (Double_t value, Double_t adjparameter) {
  value = value*gRandom->Gaus(1, adjparameter*value);
  return value;
}

 Double_t ConstFactor (Double_t value, Double_t adjparameter) {
   value = value*adjparameter;
   return value;
 }

 Double_t CurveOffset (Double_t value, Double_t adjparameter, Double_t echarge) {
   value = echarge*(1/(echarge/value + adjparameter));
   return value;
 }


//void Zanalysis (string pTparamfile="./pTparameters.csv", string rootfile="/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root") {
//void Zanalysis (vector< vector<Double_t>> pTparams, string rootfile="./13TeV_Z_PowhegPythiaDefault.root") {
void Zanalysis (vector< vector<Double_t>> pTparams, string rootfile="/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root") {

  cout << "Zanalysis has been called" << endl;

  TFile *input = TFile::Open(rootfile.c_str());
  TFile *output = new TFile("./rootfiles/Zsampletest.root", "RECREATE"); 
  TTree *MCDecayTree;
  input->GetObject("MCDecayTree", MCDecayTree);

  cout << "# of Decay tree events: " << MCDecayTree->GetEntries() << endl;

  int npTmethods = 4;
  string pTmethods[4] = {"GausSmear", "GausSmear_pTdependent", "ConstFactor", "CurveOffset"};
  //vector< vector<Double_t> > pTparams(npTmethdos);
  int hist_dims[3] = {40,80,100};

  vector< vector<TH1F *> > toys(npTmethods); 
  TH1F *nominalH = new TH1F("Znominal", "Invariant mass of muons from Z decay", hist_dims[0], hist_dims[1], hist_dims[2]);
  vector<TH1F *>::iterator toyit;
  vector<Double_t>::iterator pTparamsit;

  Double_t muMass = 0.1056583745;
  Double_t echarge = 1.602176565e-19;

  ////////////// PUT CSV READING CODE HERE /////////////////////
  /*
  ifstream input;
  input.open(pTparamfile);

  if (input.fail()) {
    cerr << "Error in pT parameter file" <<endl;
    exit(1)
  } 
  */ 
  ///////// from here assume pT params is initialized for the time being ////////////////////////////////////////////
  
  char hist_name[100];
  char hist_title[100];
  unsigned int methodindex, toyindex;

  cout << "size of pTparams: " << pTparams.size() << " size of pTparams row: " << pTparams[0].size() <<endl;
  ////////////// LOOP TO INITIALIZE THE HISTOGRAMS /////////////////////////////////////
  for (methodindex = 0; methodindex < npTmethods; ++methodindex) {
    
    for (toyindex = 0; toyindex < pTparams[methodindex].size(); ++toyindex) {
      sprintf(hist_name, "%sZ%i", pTmethods[methodindex].c_str(), toyindex);
      sprintf(hist_title, "%s: %f", pTmethods[methodindex].c_str(), pTparams[methodindex][toyindex]);
       
      cout << "Initializing histogram " << hist_name << endl;
      TH1F *hpT_toy = new TH1F(hist_name, hist_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      cout << "error marker 1" << endl;
      toys[methodindex].push_back(hpT_toy);
      cout << "error marker 2" << endl;
    }
  }

  cout << "exiting histogram initialization loop." << endl;
  //////////// toys VECTOR SHOULD NOW BE FILLED WITH EMPTY HISTOGRAMS //////////////////

  Float_t prop_M, mup_PT, mup_ETA, mup_PHI, mum_PT, mum_ETA, mum_PHI;
  Float_t mup_PTadj, mum_PTadj;

  MCDecayTree->SetBranchAddress("mup_PT", &mup_PT);
  MCDecayTree->SetBranchAddress("mum_PT", &mum_PT);
  MCDecayTree->SetBranchAddress("prop_M", &prop_M);
  MCDecayTree->SetBranchAddress("mup_ETA", &mup_ETA);
  MCDecayTree->SetBranchAddress("mup_PHI", &mup_PHI);
  MCDecayTree->SetBranchAddress("mum_ETA", &mum_ETA);
  MCDecayTree->SetBranchAddress("mum_PHI", &mum_PHI);

  cout << "error marker 3" << endl;

  MCDecayTree->SetBranchStatus("*", 0);
  MCDecayTree->SetBranchStatus("mup_PT", 1);
  MCDecayTree->SetBranchStatus("mum_PT", 1);
  MCDecayTree->SetBranchStatus("mup_ETA", 1);
  MCDecayTree->SetBranchStatus("mum_ETA", 1);
  MCDecayTree->SetBranchStatus("mup_PHI", 1);
  MCDecayTree->SetBranchStatus("mum_PHI", 1);

  //MCDecayTree->SetBranchStatus("prop_M",1);

  TLorentzVector mup, mum, musum;

  int nentries = MCDecayTree->GetEntries();
  int nbytes = 0;

  for (int eventit=0; eventit < nentries; ++eventit) {
    cout << "error marker 4" <<endl;
    nbytes += MCDecayTree->GetEntry(eventit);

    if (eventit%2 == 0){
      
      toyindex=0;
      for ( auto pTparam0 : pTparams[0] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	//cout << "error marker 4" <<endl;
	mup_PTadj = GausSmear(mup_PT, pTparam0);
	mum_PTadj = GausSmear(mum_PT, pTparam0);

	mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	musum = mup + mum;

	toys[0][toyindex]->Fill(musum.M());
	++toyindex;
      }
      
      if (npTmethods > 1) {
	toyindex=0;
	for ( auto pTparam1 : pTparams[1] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	    //cout << "error marker 4" <<endl;
	  mup_PTadj = GausSmear_pTdependent(mup_PT, pTparam1);
	  mum_PTadj = GausSmear_pTdependent(mum_PT, pTparam1);

	  mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	  mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	  musum = mup + mum;

	  toys[1][toyindex]->Fill(musum.M());
	  ++toyindex;
	}
      }
	
      if (npTmethods > 2) {
	toyindex=0;
	for ( auto pTparam2 : pTparams[2] ) {
	  mup_PTadj = ConstFactor(mup_PT, pTparam2);
	  mum_PTadj = ConstFactor(mum_PT, pTparam2);

	  mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	  mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	  musum = mup + mum;

	  toys[2][toyindex]->Fill(musum.M());
	  ++toyindex;
	}
      }
	
      if (npTmethods > 3) {
	toyindex=0;
	for ( auto pTparam3 : pTparams[3] ) {
	  mup_PTadj = CurveOffset(mup_PT, pTparam3, echarge);
	  mum_PTadj = CurveOffset(mum_PT, pTparam3, -echarge);

	  mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	  mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	  musum = mup + mum;

	  toys[3][toyindex]->Fill(musum.M());
	  ++toyindex;
	}
      }

    } else {
      mup.SetPtEtaPhiM(mup_PT, mup_ETA, mup_PHI, muMass);
      mum.SetPtEtaPhiM(mum_PT, mum_ETA, mum_PHI, muMass);

      musum = mup + mum;
      nominalH->Fill(musum.M());
    }
  }
  cout << "Histograms should all be filled" << endl;

  output->Write();
  output->Close();
}

 
int main ( ) {

  int hist_dims[3] = {40,30,50};

  int ndata_files = 1;
  string Wcharge = "Wm";
  int npTmethods = 4;

  

  int ntemplates = 11;
  int ntoys = 6;
  double pTparam_limits[8] = {0.0, 0.025, 0.0, 0.001, 0.999, 1.0005, 0.0, 0.3e-23};

  Double_t Mnom = 80.4;  
  Double_t gamma = 2.15553;
  Double_t echarge = 1.602176565e-19;
  
  

  vector<TH1F *> templates;
  vector< vector<TH1F *> > toys(npTmethods);
  vector<Double_t> Wmasses;
  vector< vector<Double_t> > pTparams(npTmethods);

  vector<TH1F *>::iterator toyit, templateit;
  vector<Double_t>::iterator Wmassit, pTparamsit;

  time_t rawtime;
  struct tm * timeinfo;
  char tbuffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(tbuffer,sizeof(tbuffer),"%d%m%Y-%I%M",timeinfo);
  string str(tbuffer);

  string fileinfo = Wcharge + "_" + to_string(ndata_files) + "ntuples_" + tbuffer;
  string output_name = "~/oxford-lhcb-wmass/rootfiles/" + fileinfo + ".root";

  TH1::SetDefaultSumw2();

  //create output file  	
  TFile *output = new TFile(output_name.c_str(),"RECREATE");

  stringstream sreweight;
  //snominal << fixed << setprecision(2) << Mnom;
  string template_name;
  string toy_name1 = "GausSmear";
  string toy_name2 = "GausSmear_pTdependent";
  string toy_name3 = "ConstFactor";
  string toy_name4 = "CurveOffset";

  cout << "error marker 1" << endl;

  // Initialize a vector of W mass hypotheses which will assist in iterating over and filling template histograms
  for (Double_t Mhyp=79.8; Mhyp<=80.8; Mhyp+=(80.8-79.8)/(ntemplates-1)) { //Mhyp for mass hypothesis
    Wmasses.push_back(Mhyp);  
    sreweight << fixed << setprecision(2) << Mhyp;
    template_name = Wcharge + sreweight.str();

    TH1F *hweighted_template = new TH1F(template_name.c_str(), "mu_PT", hist_dims[0], hist_dims[1], hist_dims[2]);    
    templates.push_back(hweighted_template);
    
    //Resetting the stringstream
    sreweight.str(string());
  }

  cout << "error marker 2" << endl;
  
  //stringstream value;
  int toyindex=0, templateindex=0;
  char hist_title[100];

  for (Double_t pTparam = pTparam_limits[0]; pTparam <= pTparam_limits[1]; pTparam += (pTparam_limits[1]-pTparam_limits[0])/(ntoys-1)) {
    pTparams[0].push_back(pTparam);
    //value << fixed << setprecision(3) << pTparam;
    toy_name1 =  toy_name1  + Wcharge + to_string(toyindex);
    sprintf(hist_title, "sigma mu_PT ~ %f", pTparam);
    TH1F *hpT_toy = new TH1F(toy_name1.c_str(), hist_title, hist_dims[0], hist_dims[1], hist_dims[2]);
    toys[0].push_back(hpT_toy);
    
    toy_name1 = "GausSmear";
    //value.str(string());
    ++toyindex;
  }

  cout << "error marker 3" << endl;

  if (npTmethods > 1) {
    toyindex = 0;
    for (Double_t pTparam = pTparam_limits[2]; pTparam <= pTparam_limits[3]; pTparam += (pTparam_limits[3]-pTparam_limits[2])/(ntoys-1)) {
      pTparams[1].push_back(pTparam);

      //value << fixed << setprecision(5) << pTparam;

      toy_name2 = toy_name2  + Wcharge + to_string(toyindex);
      sprintf(hist_title, "sigma mu_PT ~ %f * mu_PT", pTparam);
      TH1F *hpT_toy = new TH1F(toy_name2.c_str(), hist_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[1].push_back(hpT_toy);

      toy_name2 = "GausSmear_pTdependent";
      //value.str(string());
      ++toyindex;
    }
    cout << "error marker 4" << endl;
  }

  if (npTmethods > 2) {
    toyindex = 0;
    for (Double_t pTparam = pTparam_limits[4]; pTparam <= pTparam_limits[5]; pTparam += (pTparam_limits[5]-pTparam_limits[4])/(ntoys-1)) {
      pTparams[2].push_back(pTparam);

      //value << fixed << setprecision(4) << pTparam;
      toy_name3 = toy_name3 + Wcharge + to_string(toyindex);
      sprintf(hist_title, "mu_PT -> %f * mu_PT", pTparam);
      TH1F *hpT_toy = new TH1F(toy_name3.c_str(), hist_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[2].push_back(hpT_toy);

      toy_name3 = "ConstFactor";
      //value.str(string());
      ++toyindex;
    }
    cout <<"error marker 5" << endl;
  }

  if (npTmethods > 3) {
    stringstream coffset_sn;
    toyindex=0;
    for (Double_t pTparam = pTparam_limits[6]; pTparam <= pTparam_limits[7]; pTparam += (pTparam_limits[7]-pTparam_limits[6])/(ntoys-1)) {
      pTparams[3].push_back(pTparam);

      coffset_sn << scientific << setprecision(2) << pTparam;

      toy_name4 = toy_name4 + Wcharge + to_string(toyindex);
      sprintf(hist_title, "charge/mu_PT -> charge/mu_PT + %s", coffset_sn.str().c_str());
      TH1F *hpT_toy = new TH1F(toy_name4.c_str(), hist_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[3].push_back(hpT_toy);

      toy_name4 = "CurveOffset";
      coffset_sn.str(string());
      //value.str(string());
      ++toyindex;
    }
    cout << "error marker 6" << endl;
  }

  
  string parameterfile = "./pTparameters/" + fileinfo + ".csv";
  ofstream pT_ofs (parameterfile, ofstream::out);

  for (unsigned int vect_row = 0; vect_row < pTparams.size(); ++vect_row) {
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


  //  cout << "Error marker 2" << endl;

  TChain *MCDecayTree = new TChain("MCDecayTree");
  char filename[200];

  for(int k=1;k<=ndata_files;k++){
    if (k<10){
      sprintf(filename,"/data/lhcb/users/pili/forShannon/13000__%s__PowhegPythia__as0.138_IKT1.0__evts0__seed000%u.root", Wcharge.c_str(), k);
    } else sprintf(filename,"/data/lhcb/users/pili/forShannon/13000__%s__PowhegPythia__as0.138_IKT1.0__evts0__seed00%u.root", Wcharge.c_str(), k);
    cout <<"filename is  " << filename << endl;
    MCDecayTree->Add(filename);
  }
  
//Declaration of leaves types
  Float_t prop_M, mu_PT;
  Float_t mu_PTadj;
  MCDecayTree->SetBranchAddress("mu_PT", &mu_PT);
  MCDecayTree->SetBranchAddress("prop_M", &prop_M);

  MCDecayTree->SetBranchStatus("*", 0);
  MCDecayTree->SetBranchStatus("mu_PT", 1);
  MCDecayTree->SetBranchStatus("prop_M",1);

  Long64_t nentries = MCDecayTree->GetEntries();
  cout << "nentries: " << nentries << endl;
  Long64_t nbytes = 0;

  cout << "TESTING STRING COMPARE STATEMENT" << endl;
  if (Wcharge.compare("Wm") == 0) {
    cout << "The compare method correctly recognizes the Wcharge string as Wm" << endl;
  } else {
    cout << "Compare does not recognize equality of Wm with Wcharge" << endl;
  }

  TH1F *test1 = new TH1F("test1", "Fill without adjustment", hist_dims[0], hist_dims[1], hist_dims[2]);
  TH1F *test2 = new TH1F("test2", "Fill without adjustment", hist_dims[0], hist_dims[1], hist_dims[2]);

  TH1F *propM = new TH1F("propM", "propM from simulation", 40, 75, 85);
  propM->GetXaxis()->SetTitle("Mass of W propagator");
  propM->GetYaxis()->SetTitle("Counts");

  if(Wcharge.compare("Wm") == 0) {   
    echarge = -echarge;
  }

  for (Long64_t eventit=0; eventit < nentries; eventit++) { //usually i< nentries*split_ratio for full data set
    nbytes += MCDecayTree->GetEntry(eventit);
    propM->Fill(prop_M);

    if (eventit%2 == 0){
      
      test1->Fill(mu_PT);

      toyindex=0;
      for ( auto pTparam0 : pTparams[0] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	//cout << "error marker 4" <<endl;
	mu_PTadj = GausSmear(mu_PT, pTparam0);
	toys[0][toyindex]->Fill(mu_PTadj);
	++toyindex;
      }
      
      if (npTmethods > 1) {
	toyindex=0;
	for ( auto pTparam1 : pTparams[1] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	    //cout << "error marker 4" <<endl;
	  mu_PTadj = GausSmear_pTdependent(mu_PT, pTparam1);
	  toys[1][toyindex]->Fill(mu_PTadj);
	  ++toyindex;
	}
      }
	
      if (npTmethods > 2) {
	toyindex=0;
	for ( auto pTparam2 : pTparams[2] ) {
	  mu_PTadj = ConstFactor(mu_PT, pTparam2);
	  toys[2][toyindex]->Fill(mu_PTadj);
	  ++toyindex;
	}
      }
	
      if (npTmethods > 3) {
	toyindex=0;
	for ( auto pTparam3 : pTparams[3] ) {
	  
	  mu_PTadj = CurveOffset(mu_PT, pTparam3, echarge);
	  toys[3][toyindex]->Fill(mu_PTadj);
	  ++toyindex;
	}
      }

    } else {
      templateindex = 0;
      const auto denominator = TMath::BreitWigner(prop_M, Mnom, gamma); 
      
      test2->Fill(mu_PT);

      for ( auto Wmass : Wmasses ){ //Wmassit = Wmasses.begin(); Wmassit != Wmasses.end(); ++Wmassit) {
	  //cout << "error marker 5" <<endl;
	templates[templateindex]->Fill(mu_PT, TMath::BreitWigner(prop_M, Wmass, gamma)/denominator);
	  //templateindex = (templateindex + 1) % templates.size();
	++templateindex;
	  //cout << "templateindex is currently: " << templateindex <<endl;
      }
    }
  }
 

  bool firsthist = true;
  int colourit = 1;
  string templatesHname = "~/oxford-lhcb-wmass/plots/" + Wcharge + "templatesHist.png";
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
  c->Print(templatesHname.c_str());
  c->Close();    


  // PERFORMING A CHI SQUARE TEST

  cout << "All vectors and iterators initialized.  template_vect is of length " << templates.size() << ". toy_vect is of length " <<  toys.size() << ". Wmass_vect is of length " << Wmasses.size() << endl;
  cout << templates.front() << endl;
  //char scaledhname[300];

  
  TGraph *chi2Plot = new TGraph(ntemplates);
  
  char chi2plot_name[50];
  char legend_header[300];
  Double_t chi2point;
  toyindex=0;
  templateindex=0;  

  TCanvas* c2 = new TCanvas("cchi2", "chi2 with different W mass hypotheses");
  
  for (int pTmethod = 0; pTmethod < npTmethods; ++pTmethod) {
    toyindex=0;
    for (toyit = toys[pTmethod].begin(); toyit != toys[pTmethod].end(); toyit++) {

     cout << "performing loop over toys" << endl;
     if ((*toyit)->GetSumw2N() == 0) {
       cout << "Warning! Weights do not seem to be stored" << endl;
     }

     for (templateit = templates.begin(); templateit != templates.end(); templateit++){

       (*templateit)->Scale((*toyit)->Integral()/(*templateit)->Integral());

       sprintf(chi2plot_name, "chi2plot%s%u%u", Wcharge.c_str(), pTmethod, toyindex);

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
     sprintf(legend_header, "Chi2: Gaussian pT smearing %f against W mass hypotheses", pTparams[pTmethod][toyindex]);
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
  }
  c2->Close();

  
  output->Write();
  output->Close();

  Zanalysis(pTparams);

  return 0;
}

