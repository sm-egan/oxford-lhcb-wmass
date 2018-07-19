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
#include "TGraph.h"
#include <vector>

#include "TemplateStruct.h"


using namespace std;


TemplateStruct::TemplateStruct () {
  this->split_ratio = 0.5;	
  this->output_name = "~/oxford-lhcb-wmass/rootfiles/create_templates.root";
}

TemplateStruct::TemplateStruct (Double_t data_ratio, TString root_file) {
  this->split_ratio = data_ratio;
  this->output_name = root_file;
}

void TemplateStruct::drawH_BW (TChain* EventChain, string ReweightBranch, string HistBranch, int hist_dims[3], Double_t nominal_mean, Double_t reweight_mean) {
  //Set the title of the histogram to include the Wmass hypothesis
  stringstream snominal, sreweight;
  snominal << fixed << setprecision(2) << nominal_mean;
  sreweight << fixed << setprecision(2) << reweight_mean;

  string s = "Reweight"+ sreweight.str() +"Nominal"+ snominal.str();
  Double_t gamma = 2.15553;
  Long64_t maxentries = lrint((EventChain->GetEntries())*(this->split_ratio));

  TH1::SetDefaultSumw2();
  TH1F *hweighted = new TH1F(s.c_str(), HistBranch.c_str(), hist_dims[0], hist_dims[1], hist_dims[2]);
  
  //Create strings of characters to send to the draw expression which integrate the passed parameters
  char varexp[100];
  sprintf(varexp, "%s >> %s", HistBranch.c_str(), s.c_str());
  char selection[100];
  sprintf(selection, "%s > %u && %s < %u", HistBranch.c_str(), hist_dims[1], HistBranch.c_str(), hist_dims[2]);
  char weightexp[100];
  sprintf(weightexp, "(%s)*(TMath::BreitWigner(%s, %f, %f)/TMath::BreitWigner(%s,%f,%f))", 
	  selection,ReweightBranch.c_str(),reweight_mean,gamma,ReweightBranch.c_str(),nominal_mean,gamma);

  //  EventChain->Draw(varexp, weightexp,"", 1000); 
  EventChain->Draw(varexp, weightexp,"", maxentries, maxentries); 
  this->templates.push_back(hweighted);
  this->Wmasses.push_back(reweight_mean);
}

void TemplateStruct::create_templates(int hist_dims[3], int ntemplates=5, int ndata_files=2) {

  TH1::SetDefaultSumw2();

  //define histogram
  TH1F *h_muPT= new TH1F("h_mu_PT","#mu P_{T}",hist_dims[0], hist_dims[1], hist_dims[2]);

  //create output file  	
  TFile *output = new TFile(this->output_name,"RECREATE");
 
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
   // ask Martina what this is for - I don't see it being used
   Long64_t nbytes = 0;

   for (Long64_t i=0; i<(nentries*split_ratio);i++) { //usually i< nentries*split_ratio for full data set

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
      
      /* The goal of the lines commented out here was to make sure that there would always be n templates, even when the increment of the addition does not divide perfectly into 1
      if (((Mhyp +(80.8-79.8)/(ntemplates-1)) > 80.8) && (Mhyp != 80.8)) {
	this->drawH_BW(MCDecayTree, "prop_M", "mu_PT", hist_dims, Mnom, 80.8);
	}*/
    }

    
    //TemplateStruct ts; This function is already acting on an instance
    //ts.Wmasses = Wmass_vect; Should be filled already in the drawH_BW function
    //ts.templates = hist_vect; "      
    //ts.toys.push_back(h_muPT); 

    // Add the single toy histogram to the TemplateStruct vector.  Eventually this will be modified to create a loop
    this->toys.push_back(h_muPT);
    //output->Close();

    //this->template_chi2();
    
    vector<TH1F *>::iterator toyit;

  cout << "All vectors and iterators initialized.  template_vect is of length " << (this->templates).size() << ". toy_vect is of length " <<  (this->toys).size() << ". Wmass_vect is of length " << (this->Wmasses).size() << endl;
  cout << (this->templates).front() << endl;
  char scaledhname[200];

  int toyindex=0, templateindex=0;
  // Declare a vector of chi-square results, where the first dim represents the toys and second dim gives result for each template
  //vector<Double_t> chi2_results(template_vect.size());
  //Double_t chi2_results[toy_vect.size()][template_vect.size()];
  //  Double_t chi2_results[1][5];
  
  // Copy template hypothesis info onto an array of fixed size to hopefully be accepted by TGraph->DrawGraph()
  
  Double_t Wmass_arr[ntemplates];
  cout << "Copying TemplateStruct->Wmasses to an array" << endl;
  copy((this->Wmasses).begin(), (this->Wmasses).end(), Wmass_arr);
  
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

  for (toyit = (this->toys).begin(); toyit != (this->toys).end(); toyit++) {
    //(*toyit)->Scale(1 / ((*toyit)->Integral()), "width");
    //sprintf(scaledhname, "%s_scaled", (*toyit)->GetName());
    //output->WriteTObject(*toyit, scaledhname, "Overwrite");
    cout << "performing loop over toys" << endl;
    if ((*toyit)->GetSumw2N() == 0) {
      cout << "Warning! Weights do not seem to be stored" << endl;
    }

    for (vector<TH1F *>::iterator templateit = (this->templates).begin(); templateit != (this->templates).end(); templateit++){
      (*toyit)->Scale((*templateit)->Integral()/(*toyit)->Integral());
      sprintf(scaledhname, "%s_scaledto_%s", (*toyit)->GetName(), (*templateit)->GetName());
      output->WriteTObject(*toyit, scaledhname, "Overwrite");
      chi2point = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
      
      cout << "Copying chi2 result for toy " << toyindex << " and template " << templateindex << ": " << chi2point << endl;
      chi2Plot->SetPoint(templateindex,Wmass_arr[templateindex],chi2point);
      //chi2_results[toyindex][templateindex] = (*templateit)->Chi2Test((*toyit), "WW");
      //cout << "Copying dummy result to chi2_results array for testing" << endl;
      //chi2_results[toyindex][templateindex] = 12.0;
      // Increment the template index, but adjust with modulo to avoid segmentation fault 
      templateindex = (templateindex + 1) % (this->templates).size();
    }
    cout << "Drawing graph of toy " << toyindex << endl;
    //chi2Plot->DrawGraph((this->templates).size(), Wmass_arr, chi2_results[toyindex]);
    chi2Plot->SetMarkerColor(4);
    chi2Plot->SetMarkerStyle(8);
    chi2Plot->Draw("AP");
    //output->WriteTObject(TGraph(this->Wmasses, chi2_results), 
    //graph_name.c_str(), "Overwrite");
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

  output->Close();
}

void TemplateStruct::plot_overlaid_histograms () {
  TFile *output = TFile::Open(this->output_name, "UPDATE");  
  bool firsthist = true;
  int colourit = 1;
  vector<TH1F *>::iterator templateit;
    // Below might be necessary for larger number of hypotheses
    //int lineit = 1;
    
  TCanvas* c = new TCanvas("cBW", "mu PT with different W mass hypotheses");
      
  for (templateit = this->templates.begin(); templateit != this->templates.end(); templateit++, colourit++) {
      //brackets around *templateit ensure that we are acting on the pointer to the TH1F, not the iterative pointer to the pointer  
    output->WriteTObject(*templateit, (*templateit)->GetName(),"Overwrite");
    (*templateit)->SetLineColor(colourit);

    if(firsthist) {
      (*templateit)->Draw();
      firsthist = false;
    } else (*templateit)->Draw("SAME");
  }
  c->Print("~/oxford-lhcb-wmass/plots/templatesHist.png");
  c->Print("~/oxford-lhcb-wmass/plots/templatesHist.pdf");
  c->Close();    
}

void TemplateStruct::template_chi2 () {
  
  TFile *output = TFile::Open(this->output_name, "UPDATE");  
  vector<TH1F *>::iterator templateit;
  vector<TH1F *>::iterator toyit;
  int ntemplates = this->Wmasses.size();

  cout << "All vectors and iterators initialized.  template_vect is of length " << (this->templates).size() << ". toy_vect is of length " <<  (this->toys).size() << ". Wmass_vect is of length " << (this->Wmasses).size() << endl;
  cout << (this->templates).front() << endl;
  char scaledhname[200];

  int toyindex=0, templateindex=0;  
  // Copy template hypothesis info onto an array of fixed size to hopefully be accepted by TGraph->DrawGraph()
  
  Double_t Wmass_arr[ntemplates];
  cout << "Copying TemplateStruct->Wmasses to an array" << endl;
  copy((this->Wmasses).begin(), (this->Wmasses).end(), Wmass_arr);
  
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

  for (toyit = (this->toys).begin(); toyit != (this->toys).end(); toyit++) {
    //(*toyit)->Scale(1 / ((*toyit)->Integral()), "width");
    //sprintf(scaledhname, "%s_scaled", (*toyit)->GetName());
    //output->WriteTObject(*toyit, scaledhname, "Overwrite");
    cout << "performing loop over toys" << endl;
    if ((*toyit)->GetSumw2N() == 0) {
      cout << "Warning! Weights do not seem to be stored" << endl;
    }

    for (templateit = (this->templates).begin(); templateit != (this->templates).end(); templateit++){
      (*toyit)->Scale((*templateit)->Integral()/(*toyit)->Integral());
      sprintf(scaledhname, "%s_scaledto_%s", (*toyit)->GetName(), (*templateit)->GetName());
      output->WriteTObject(*toyit, scaledhname, "Overwrite");
      chi2point = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
      
      cout << "Copying chi2 result for toy " << toyindex << " and template " << templateindex << ": " << chi2point << endl;
      chi2Plot->SetPoint(templateindex,Wmass_arr[templateindex],chi2point);
      //chi2_results[toyindex][templateindex] = (*templateit)->Chi2Test((*toyit), "WW");
      //cout << "Copying dummy result to chi2_results array for testing" << endl;
      //chi2_results[toyindex][templateindex] = 12.0;
      // Increment the template index, but adjust with modulo to avoid segmentation fault 
      templateindex = (templateindex + 1) % (this->templates).size();
    }
    cout << "Drawing graph of toy " << toyindex << endl;
    //chi2Plot->DrawGraph((this->templates).size(), Wmass_arr, chi2_results[toyindex]);
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

  output->Close();

  /*  cout << "All vectors and iterators initialized.  template_vect is of length " << (this->templates).size() << ". toy_vect is of length " <<  (this->toys).size() << ". Wmass_vect is of length " << (this->Wmasses).size() << endl;
  cout << (this->templates).front() << endl;
  
  //Unit normalize all of the template histograms

  for (templateit = (this->templates).begin(); templateit != (this->templates).end(); templateit++) {
    cout << "Scaling template histogram at index " << distance((this->templates).begin(),templateit) << endl;
    (*templateit)->Scale(1/((*templateit)->Integral()));
    cout << "histogram has been scaled." << endl;
  }
  cout << "Template histograms have been scaled." << endl;

  int toyindex=0, templateindex=0;
  // Declare a vector of chi-square results, where the first dim represents the toys and second dim gives result for each template
  //vector<Double_t> chi2_results(template_vect.size());
  //Double_t chi2_results[toy_vect.size()][template_vect.size()];
  Double_t chi2_results[1][5];
  
  // Copy template hypothesis info onto an array of fixed size to hopefully be accepted by TGraph->DrawGraph()
  Double_t Wmass_arr[5];
  cout << "Copying TemplateStruct->Wmasses to an array" << endl;
  copy((this->Wmasses).begin(), (this->Wmasses).end(), Wmass_arr);
  TGraph *chi2Plot = new TGraph();
  string graph_name = "uniquename";

//Unit normalize all of the toy vectors and perform chi square test, saving the data to an array as you go
  TCanvas* c2 = new TCanvas("cchi2", "chi2 with different W mass hypotheses");

  for (toyit = (this->toys).begin(); toyit != (this->toys).end(); toyit++) {
    (*toyit)->Scale(1 / ((*toyit)->Integral()));
    for (templateit = (this->templates).begin(); templateit != (this->templates).end(); templateit++){
      
      cout << "Copying chi2 result for toy " << toyindex << " and template " << templateindex << endl;
      chi2_results[toyindex][templateindex] = (*templateit)->Chi2Test((*toyit), "WW");
      //cout << "Copying dummy result to chi2_results array for testing" << endl;
      //chi2_results[toyindex][templateindex] = 12.0;
      // Increment the template index, but adjust with modulo to avoid segmentation fault 
      templateindex = (templateindex + 1) % (this->templates).size();
    }
    cout << "Drawing graph of toy " << toyindex << endl;
    chi2Plot->DrawGraph((this->templates).size(), Wmass_arr, chi2_results[toyindex]);
    
    //output->WriteTObject(TGraph(this->Wmasses, chi2_results), 
    //graph_name.c_str(), "Overwrite");
    cout << "Writing chi2 plot to root file." << endl;
    output->WriteTObject(chi2Plot, graph_name.c_str(), "Overwrite");
    toyindex++;
  }
  cout << "Saving chi2 plot as an image." << endl;
  c2->Print("./plots/chi2Plot.png");
  c2->Close(); */ 
        // A potential next step: open TCanvas to write TGraphs of different toys to the same plot

}

/*
int main () {
  return 0;
}
*/
