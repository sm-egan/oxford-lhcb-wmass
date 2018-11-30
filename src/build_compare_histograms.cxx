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
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include "TSystem.h"
#include <vector>
#include "TemplateStruct.h"
#include "TRandom.h"
#include <fstream>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include <TLorentzVector.h>
#include <cmath>

///data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root
using namespace std;
namespace po = boost::program_options;

Double_t GausSmear (Double_t pT, Double_t adjparameter, bool generate_smear=true) {
  if (generate_smear) {
    pT = pT*gRandom->Gaus(1, adjparameter);
  } else {
    //cout << "GausSmear: Multiplying by pT smear parameter of " << adjparameter << endl;
    pT = pT*adjparameter;
  }
  return pT;
}

Double_t GausSmear_pTdependent (Double_t pT, Double_t adjparameter, bool generate_smear=true) {
  if (generate_smear) {
    pT = pT*gRandom->Gaus(1, adjparameter*pT);
  } else {
    pT = pT*adjparameter;
  }
  return pT;
}

Double_t GausSmear_pdependent (Double_t pT, Double_t adjparameter, Double_t eta, bool generate_smear=true) {
  if (generate_smear) {
    pT = pT*gRandom->Gaus(1, adjparameter*pT*cosh(eta));
  } else {
    //cout << "GausSmear_pdependent: Multiplying by pT smear parameter of " << adjparameter << endl;
    pT = pT*adjparameter;
  }
  return pT;
}

Double_t MomentumScale (Double_t pT, Double_t adjparameter) {
  pT = pT*(adjparameter);
  return pT;
}

Double_t CurvatureBias (Double_t pT, Double_t adjparameter, Double_t echarge) { //, Double_t eta) {
  pT = (echarge/(echarge/(pT) + adjparameter));
  return pT;
}

Double_t CurvatureBias (Double_t pT, Double_t adjparameter, Double_t echarge, Double_t eta) {
  pT = (echarge/(echarge/(pT*cosh(eta)) + adjparameter))/cosh(eta);
  return pT;
}

Double_t ResolutionSmear (Double_t pT, Double_t constcoeff, Double_t lincoeff, Double_t eta, bool generate_smear=true) {
  pT = GausSmear(pT, constcoeff, generate_smear);
  pT = GausSmear_pdependent(pT, lincoeff, eta, generate_smear);
  return pT;
}

string make_timestamp (){
  time_t rawtime;
  struct tm * timeinfo;
  char tbuffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(tbuffer,sizeof(tbuffer),"%d%m%Y-%H%M",timeinfo);
  string str(tbuffer);
  return tbuffer;
}

string trunc_double (double x, int precision, string format = "decimal") {
  stringstream ss;
  //cout << "trunc_double input: " << x << endl;
  if (format.compare("decimal") == 0) {
    ss << fixed << setprecision(precision) << x;
  } else if (format.compare("scientific") == 0) {
    ss << scientific << setprecision(precision) << x;
  }
  string outstr = ss.str();
  //cout << "trunc_double stringstream: " << outstr << endl;
  return outstr;
}
/*
//Doing this as a helper method probably won't work because histograms are out of scope and we won't have access to any of the histogram names in order to check which x section to use

Double_t Chi2Test_lumiscaling(TH1F templateH, TH1F toyH) {

  Double_t lhcb_luminosity = 6.0; //Predicted fb^-1 at the end of run 2
  Double_t xs_Zmumu = 198000.0, xs_Wpmunu = 1093600.0, xs_Wmmunu = 818400;

}
*/

//void Zanalysis (vector< vector<Double_t>> pTparams, string rootfile="./13TeV_Z_PowhegPythiaDefault.root") {

vector<Double_t> scale_vector (vector<Double_t> to_scale, Double_t factor) {
  for (vector<Double_t>::iterator scaleit = to_scale.begin(); scaleit != to_scale.end(); ++scaleit) {
    cout << "Scaling value " << (*scaleit);
    (*scaleit) *= factor;
    cout << "to " << (*scaleit) << endl;
  }
  return to_scale;
}

template <typename TData>
void tcanvas_from_vector (vector<TData> vector_to_plot, int step = 1, bool add_legend = false, const vector<Double_t> legend_info = vector<Double_t>(), bool add_title = true, string canvas_name="", string canvas_title="", string plot_name="") {
  bool firsthist = true;
  int colourit = 0;
  string plot_namepng;
  // Flag which forces return if you give only one of the strings.  The flag works off the assumption that you have to upate profile_title if you want to update plot_name, so it tests only if plot_name is still
  /*
  if ((canvas_name.compare("") != 0) && (plot_name.compare("") == 0)) {
    cout << "tcanvas_from_vector() requires that either all string fields are specified, or that none of them are.  Please revise passage of arguments" << endl;
    return;
  }
  */
 
  if (canvas_name.compare("") == 0) {
    string c_name((*vector_to_plot.begin())->GetName());
    canvas_name = c_name.substr(0, c_name.length() -1) + "_combined";
  }
  
  if (canvas_title.compare("") == 0) {
    string c_title((*vector_to_plot.begin())->GetTitle());
    canvas_title = c_title;
  } 

  if (plot_name.compare("") == 0) {
   plot_name = "/home/egan/oxford-lhcb-wmass/plots/" + canvas_name + ".pdf";
   plot_namepng = "/home/egan/oxford-lhcb-wmass/plots/" + canvas_name + ".png";
  }
   
  int vectindex = 0;

  TCanvas *c = new TCanvas(canvas_name.c_str()); 
  if (add_title) {
    c = new TCanvas(canvas_name.c_str(), canvas_title.c_str()); 
  }

  int colours[8] = {600, 807, 1, 616, 843, 632, 434, 402};
  //gStyle->SetPalette(7, colours);
  /*
kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
  */

  if (!add_legend) {
    for (typename vector<TData>::iterator vectorit = vector_to_plot.begin(); vectorit != vector_to_plot.end(); ++vectorit) {
      if (vectindex % step == 0) {
	//(*vectorit)->SetLineColor();
	
	if(firsthist) {
	  (*vectorit)->Draw();
	  firsthist = false;
	} else (*vectorit)->Draw("SAME");
	++colourit;
      }
      ++vectindex;
    }
  } else {
    TLegend *legend = new TLegend(0.3, 0.1, 0.7, 0.5);
    legend->SetHeader("Curvature Bias Parameter");
    char leglabel[100];

    /*
    for(int vindex; vindex < vector_to_plot.size(); ++vindex) {
	if (vindex % step == 0) {
	  vector_to_plot[vindex]->SetLineColor(colourit);
	  legend->AddEntry(vector_to_plot[vindex], trunc_double(legend_info[vindex], 3, "scientific").c_str(), "l");

	  if(firsthist) {
	    vector_to_plot[vindex]->Draw();
	    firsthist = false;
	  } else vector_to_plot[vindex]->Draw("SAME");
	  ++colourit;
	}
    }
    */
    
    for (typename vector<TData>::iterator vectorit = vector_to_plot.begin(); 
	 vectorit != vector_to_plot.end(); ++vectorit) {
      if (vectindex % step == 1) {
	(*vectorit)->SetLineColor(colours[colourit]);
	legend->AddEntry((*vectorit), trunc_double(legend_info[vectindex], 3, "scientific").c_str(), "l");

	if(firsthist) {
	  (*vectorit)->Draw();
	  firsthist = false;
	} else (*vectorit)->Draw("SAME");
	++colourit;
      }
      ++vectindex;
    }
    legend->Draw();
  }

  c->Print(plot_name.c_str());
  c->Print(plot_namepng.c_str());
  c->Close();
}

void draw_colourhist2d (TH2F *hist, bool logx=false, bool logy=false, string plot_path = "/home/egan/oxford-lhcb-wmass/plots/kinematics/") {
  string plot_name(hist->GetName());   
  // The order of the following two lines is important - otherwise you will wrap the png filename around the pdf filename
  string plot_namepng = plot_path + plot_name + ".png";
  plot_name = plot_path + plot_name + ".pdf";
  //Double_t w = 500;
  //Double_t h = 500;

  TCanvas *c = new TCanvas("c", hist->GetTitle());
  gStyle->SetPalette(57);
  gStyle->SetOptStat(0);
  hist->SetContour(1000);
  
  if (logx) {
    c->SetLogx();
  }
  if (logy) {
    c->SetLogy();
  }

  hist->Draw("colz");
  //c->SetWindowSize(1.1*c->GetWw(), c->GetWh());
  c->Print(plot_name.c_str());
  c->Print(plot_namepng.c_str());
  c->Close();
}


void dimuon_analysis (vector< vector<Double_t>> pTparams, 
		      //string output_name, //= "./rootfiles/Zsampletest_lhcbcuts.root", 
		      string dimu_propagator,
		      string rootfile="",
		      bool smear_template = true,
		      bool fill_curveobs = false,
		      bool fill_profiles = false,
		      bool fill_asymsubsets = false,
		      bool fill_mscbtoys = false,
		      bool use_all_events = false) { //="/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root") {

  cout << "dimuon_analysis has been called" << endl;

  int npTmethods = 4;
  string pTmethods[4] = {"GausSmear", "GausSmear_pTdependent", "MomentumScale", "CurvatureBias"};
  int ncurveobs = 2;
  string curveobservables[2] = {"pasym", "deltamuP"};// "asym_dmuPT", "asym_dmuP"};
  int nprofiles = 2;
  string profilevars[2] = {"pasym", "deltamuP"};
    //vector< vector<Double_t> > pTparams(npTmethdos);
  string propagator;
  double hist_lims[5][2];
  int nbins = 40; 
  
  double min_mupT;
  
  if (dimu_propagator.compare("Z") == 0) { 
    cout << "Z boson sample found" << endl; 
    //bool isZ = true;
    propagator = "Z";    
    min_mupT = 20;

    hist_lims[0][0] = 80;
    hist_lims[0][1] = 100;

  } else if (dimu_propagator.find("Y") == 0) {
    
    cout << "Upsilon sample found" << endl;
    propagator = "Upsilon";
    min_mupT = 2;
    //int upsilondims[3] = {9,11,50};
    //copy(upsilondims, upsilondims+3, hist_dims.begin());
    //hist_dims = set_hist_dims(hist_dims, 9,11,25);
    nbins = 70;
    // set hist_lims[0] >= hist_lims[1] to get automatic determination of the limits
    hist_lims[0][0] = 9.35;
    hist_lims[0][1] = 9.55; 

    hist_lims[1][0] = -1;
    hist_lims[1][1] = 1;

    hist_lims[2][0] = -200;
    hist_lims[2][1] = 200;

    hist_lims[3][0] = 1;
    hist_lims[3][1] = 0;

    hist_lims[4][0] = 1;
    hist_lims[4][1] = 0;    

    string parameterfile = "./pTparameters/" + propagator + "sample.csv";
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
  
  } else {
    int hist_dims[3] = {0,100,100};
    cout << "WARNING: could not identify particle type by rootfile string" << endl;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// PREPARE TTREE OR TCHAIN OF EVENTS ///////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  TChain *MCDecayTree = new TChain("MCDecayTree");
  char filename[200];
  if (rootfile.compare("") == 0) {
    
    if (dimu_propagator.compare("Y") == 0) {
      cout << "Loading default Upsilon sample root files" << endl;
      
      for(int k=1; k<=5; k++){
	if (k<10){
	  sprintf(filename,"/home/egan/oxford-lhcb-wmass/Ysamples/default/pythia_upsilon1S_muonly_10000000events_seed00%u.root", k);
	} else sprintf(filename,"/home/egan/oxford-lhcb-wmass/Ysamples/default/pythia_upsilon1S_muonly_10000000events_seed0%u.root", k);
	cout <<"filename is  " << filename << endl;
	MCDecayTree->Add(filename);
      }

    } else if (dimu_propagator.compare("Z") == 0) {
      cout << "Loading default Z sample root file" << endl;
      //TFile *input = TFile::Open("/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root");
      //input->GetObject("MCDecayTree", MCDecayTree);
      sprintf(filename, "/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root");
      MCDecayTree->Add(filename);
    }
    //if (dimu_propagator.compare(""))
  } else {
    MCDecayTree->Add(rootfile.c_str());
  }

  cout << "# OF DECAY TREE EVENTS: " << MCDecayTree->GetEntries() << endl;

  char output_name[100];
  string tbuffer = make_timestamp();
  sprintf(output_name, "./rootfiles/%ssample_%s.root", propagator.c_str(), tbuffer.c_str());
  TFile *output = new TFile(output_name, "RECREATE"); 
 
  TH1::SetDefaultSumw2();
  char nominalH_name[50];
  sprintf(nominalH_name, "%stemplate", propagator.c_str());
  TH1F *nominalH = new TH1F(nominalH_name, "Invariant mass of muons", nbins, hist_lims[0][0], hist_lims[0][1]);

  vector< vector<TH1F *> > toys(npTmethods);
  vector< vector <TH1F *> > curveobs(2);
  //vector<vector<TProfile *> > profiles(2);
  vector <TProfile *> profiles;
  vector< vector <TH1F *> > asym_subsets(11);
    // Combined momentum scale and curvature bias toys
  vector< vector < vector <TH1F *> > > mscbtoys(11, vector<vector<TH1F *>>(11));
  
  vector<TH1F *>::iterator toyit;
  vector<Double_t>::iterator pTparamsit;

  Double_t muMass = 0.1056583745;
  Double_t echarge = 1.602176565e-19;
  //Double_t lhcb_luminosity = 6.0; //Predicted fb^-1 at the end of run 2  
  //Double_t xs_Zmumu = 198000.0; // Units in fb

  char hist_name[100];
  char cohist_name[100];
  char hist_title[100];
  char profile_name[100];
  char asymsubset_name[100];
  char asymsubset_title[100];
  unsigned int methodindex, toyindex, coindex, profindex, asymindex, cbindex;

  cout << "size of pTparams: " << pTparams.size() << " size of pTparams row: " << pTparams[0].size() <<endl;
  cout << "size of curveobs: " << curveobs.size() << endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// LOOP TO INITIALIZE THE HISTOGRAMS /////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  for (methodindex = 0; methodindex < pTparams.size(); ++methodindex) {
    for (toyindex = 0; toyindex < pTparams[methodindex].size(); ++toyindex) {
      sprintf(hist_name, "%s%s%i", pTmethods[methodindex].c_str(), propagator.c_str(), toyindex);
      sprintf(hist_title, "%s: %f", pTmethods[methodindex].c_str(), pTparams[methodindex][toyindex]);

      if ((methodindex == 2) && fill_asymsubsets && (!fill_mscbtoys)) {
	asymindex = 0;
	for (vector<vector<TH1F *>>::iterator asymit = asym_subsets.begin(); asymit != next(asym_subsets.begin(), 3); ++asymit) {
          sprintf(asymsubset_name, "%s%s_asymubset%u%u", pTmethods[methodindex].c_str(), propagator.c_str(), asymindex, toyindex);
	  sprintf(asymsubset_title, "Dimuon invariant mass for p asymmetry subset %u", asymindex);
	  TH1F *asymhist = new TH1F(asymsubset_name, asymsubset_title, 50, hist_lims[0][0], hist_lims[0][1]);
	  (*asymit).push_back(asymhist);
	  ++asymindex;
	}
      }	
      	
      if ((methodindex == 3) && (fill_curveobs || fill_profiles || fill_asymsubsets || fill_mscbtoys)) {
	sprintf(hist_title, "%s: %s", pTmethods[methodindex].c_str(), trunc_double(pTparams[methodindex][toyindex], 2, "scientific").c_str());

	if (fill_curveobs) {
	  coindex = 0;
	  for (vector< vector<TH1F *>>::iterator coit = curveobs.begin(); coit != curveobs.end(); ++coit) {	    
	    sprintf(cohist_name, "%s_%s", curveobservables[coindex].c_str(), hist_name);
	    TH1F *cohist = new TH1F(cohist_name, hist_title, nbins, hist_lims[coindex+1][0], hist_lims[coindex+1][1]);	  
	    //TH1F *cohist = new TH1F(cohist_name, hist_title, nbins, 1, 0);	  
	    (*coit).push_back(cohist);
	    ++coindex;
	  }
	}

	if (fill_profiles) {
	  profindex = 0;
	  // Require odd numbers so that just 5 curves are plotted - 11 too crowded 
	  //if ((toyindex % 2) == 1) {
	  sprintf(profile_name, "%s_%s_profile%u", propagator.c_str(), profilevars[profindex].c_str(), toyindex);
	  TProfile *prof = new TProfile(profile_name, "Invariant mass profile for dimuon system", 30, hist_lims[1][0], hist_lims[1][1],9.4,9.5);
          prof->SetStats(kFALSE);
          prof->GetXaxis()->SetTitle("Total momentum asymmetry");
          prof->GetYaxis()->SetTitle("Mean invariant mass (GeV)");
          profiles.push_back(prof);
	    //}
	  /* These lines of code should be used when profiles is a 2D vector such that both the p asymmetry and p difference are plotted
	  for (vector< vector<TProfile *>>::iterator profit = profiles.begin(); profit != profiles.end(); ++profit) {	    
	    sprintf(profile_name, "%s_%s_profile%u", propagator.c_str(), profilevars[profindex].c_str(), toyindex);
	    if (profilevars[profindex].compare("pasym") == 0) {
	      TProfile *prof = new TProfile(profile_name, "Profile of Dimuon invariant mass vs mu p asymmetry", 30, hist_lims[1][0], hist_lims[1][1]);
	      (*profit).push_back(prof);
	    } else {
	      TProfile *prof = new TProfile(profile_name, "Profile of Dimuon invariant mass vs difference in mu p", 30, hist_lims[2][0], hist_lims[2][1]);
	      (*profit).push_back(prof);
	    }
	    ++profindex;
	  */  
	}

	if (fill_asymsubsets && (!fill_mscbtoys)) {
	  asymindex=3;
	  for (vector<vector<TH1F *>>::iterator asymit = next(asym_subsets.begin(), 3); asymit != asym_subsets.end(); ++asymit) {
	    sprintf(asymsubset_name, "%s%s_asymubset%u%u", pTmethods[methodindex].c_str(), propagator.c_str(), asymindex, toyindex);
	    sprintf(asymsubset_title, "Dimuon invariant mass for p asymmetry subset %u", asymindex);
	    TH1F *asymhist2 = new TH1F(asymsubset_name, asymsubset_title, 50, hist_lims[0][0], hist_lims[0][1]);
	    (*asymit).push_back(asymhist2);
	    ++asymindex;
	  }
	} else if (fill_mscbtoys) {
	  // We have the overarching toy index that we can use to access the first dimension of mscbtoys. The goal is that after the following loop we should be able to push a 2D
	  // 11 x 3 vector onto the mcbtoys[toyindex]
	  cbindex = 0;
	  for (vector<vector<TH1F *>>::iterator cbasymit = mscbtoys[toyindex].begin(); cbasymit != mscbtoys[toyindex].end(); ++cbasymit) {
	    //cout << "The type of *cbasymit is " << typeid(*cbasymit).name() << endl;
	    for (int cbasymindex = 0; cbasymindex < 3; ++cbasymindex) {
	      sprintf(asymsubset_name, "mscb%s_asymsubset%u%u%u", propagator.c_str(), toyindex, cbindex, cbasymindex);
	      sprintf(asymsubset_title, "Momentum scale: %f, Curvature bias: %f, Asym subset: %u", pTparams[2][toyindex], pTparams[3][cbindex], cbasymindex);
	      cout << "Initializing histogram " << asymsubset_name << endl;
	      TH1F *cbasymhist = new TH1F(asymsubset_name, asymsubset_title, nbins, hist_lims[0][0], hist_lims[0][1]);
	      (*cbasymit).push_back(cbasymhist);
	      cout << asymsubset_name << " successfully added to vector" << endl;
	    }
	    mscbtoys[toyindex][cbindex] = (*cbasymit);
	    ++cbindex;
	  }	  
	}
      }
      cout << "Initializing histogram " << hist_name << endl;
      TH1F *hpT_toy = new TH1F(hist_name, hist_title, nbins, hist_lims[0][0], hist_lims[0][1]);
      toys[methodindex].push_back(hpT_toy);
    }
  }	   

  /*
  if (fill_mscbtoys) {	  
    // LOOP OVER THE MOMENTUM SCALE
    msindex = 0;
    for (vector<vector<vector<TH1F *>>>::iterator msit = mscbtoys.begin(); msit != mscbtoys.end(); ++msit) {
      cbindex = 0;
      for (vector<vector<TH1F *>>::iterator cbit = (*msit).begin(); cbit != (*msit).end(); ++cbit) {
	sprintf(mscbtoy_name, "mscb%s_asymsubset%u%u%u");
      }
      ++msindex;
    }
  }
  */

	   
  cout << "exiting histogram initialization loop." << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// FILLING ALL HISTOGRAMS AND PROFILES ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Float_t mup_PT, mup_ETA, mup_PHI, mum_PT, mum_ETA, mum_PHI;
  Float_t mup_PTadj, mum_PTadj, mup_PTsmear, mum_PTsmear;

  MCDecayTree->SetBranchAddress("mup_PT", &mup_PT);
  MCDecayTree->SetBranchAddress("mum_PT", &mum_PT);
  MCDecayTree->SetBranchAddress("mup_ETA", &mup_ETA);
  MCDecayTree->SetBranchAddress("mup_PHI", &mup_PHI);
  MCDecayTree->SetBranchAddress("mum_ETA", &mum_ETA);
  MCDecayTree->SetBranchAddress("mum_PHI", &mum_PHI);

  MCDecayTree->SetBranchStatus("*", 0);
  MCDecayTree->SetBranchStatus("mup_PT", 1);
  MCDecayTree->SetBranchStatus("mum_PT", 1);
  MCDecayTree->SetBranchStatus("mup_ETA", 1);
  MCDecayTree->SetBranchStatus("mum_ETA", 1);
  MCDecayTree->SetBranchStatus("mup_PHI", 1);
  MCDecayTree->SetBranchStatus("mum_PHI", 1);

  //////// Initialize 2D kinematic histogram /////////////
  // Notice that the current setting gives automatic ranges - maybe change this to include only LHCb acceptances
  //string dimumass2d_name = propagator + "kin2d_dimumass_eta";
  //TH2F *dimumass2d = new TH2F(dimumass2d_name.c_str(), "Histogram of dimuon system invariant mass and rapidity", 10, 1, 0, 10, 1, 0);
  //dimumass2d->GetXaxis()->SetTitle("Invariant mass (GeV)");
  //dimumass2d->GetYaxis()->SetTitle("Rapidity");

  string kinematics2d_name = propagator + "kin2d_logpT_eta";
  TH2F *kinematics2d = new TH2F(kinematics2d_name.c_str(), "Histogram of dimuon system pT and rapidity", 10, 1, 0, 10, 2, 5);
  kinematics2d->GetXaxis()->SetTitle("log(pT) (GeV)");
  kinematics2d->GetYaxis()->SetTitle("Rapidity");

  TLorentzVector mup, mum, musum;
  int nentries = MCDecayTree->GetEntries();
  int nbytes = 0;
  int mscbindex;
  Double_t asym_limit = 1.0/3.0;
  Double_t mup_constsmear, mum_constsmear, mup_psmear, mum_psmear;


  // Initializing event loop to fill histograms
  for (int eventit=0; eventit < nentries; ++eventit) {
    nbytes += MCDecayTree->GetEntry(eventit);

    // The terms from the pT parameter matrix are used to set the intermediate smear value.  The last entry in the row
    // is the highest possible smear value in range that begins at 0, so half this value is the midpoint
    mup_constsmear = gRandom->Gaus(1, 0.5*pTparams[0].back());
    mum_constsmear = gRandom->Gaus(1, 0.5*pTparams[0].back());
    mup_psmear = gRandom->Gaus(1, 0.5*pTparams[1].back()*mup_PT*cosh(mup_ETA));
    mum_psmear = gRandom->Gaus(1, 0.5*pTparams[1].back()*mum_PT*cosh(mum_ETA));

    if (smear_template) {
      mup_PTsmear = ResolutionSmear(mup_PT, mup_constsmear, mup_psmear, mup_ETA, false);
      mum_PTsmear = ResolutionSmear(mum_PT, mum_constsmear, mum_psmear, mum_ETA, false);
    } else {
      mup_PTsmear = mup_PT;
      mum_PTsmear = mum_PT;
    }

    // Only accept muons that are travelling in same direction
    if (mup_ETA*mum_ETA > 0) {
      mup_ETA = abs(mup_ETA);
      mum_ETA = abs(mum_ETA);
      // min_muPT is set according to the propagator used (upsilon or W).  This statement skips the event if the minimum Pt requirement
      // is not met or if any of the muons are outside of LHCb's acceptance range
      if ((mup_PTsmear < min_mupT) || (mum_PTsmear < min_mupT) || (mup_ETA < 2.0) || (mup_ETA > 4.5) || (mup_ETA < 2.0) || (mum_ETA > 4.5)) {
	continue;
      }
    } else continue;

    // Plot a 2D histogram of the kinematic variables WITHOUT smearing
    mup.SetPtEtaPhiM(mup_PT, mup_ETA, mup_PHI, muMass);
    mum.SetPtEtaPhiM(mum_PT, mum_ETA, mum_PHI, muMass);
    musum = mup + mum;

    //dimumass2d->Fill(musum.M(), musum.Rapidity());
    kinematics2d->Fill(log10(musum.Pt()), musum.Rapidity());

    if ((eventit%2 == 0) || use_all_events) {
      
      toyindex=0;
      
      for ( auto pTparam0 : pTparams[0] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	mup_PTadj = GausSmear(GausSmear_pdependent(mup_PT, mup_psmear, mup_ETA, false), pTparam0);
	//mup_PTadj = GausSmear(mup_PT, *pTparamsit);
	mum_PTadj = GausSmear(GausSmear_pdependent(mum_PT, mum_psmear, mum_ETA, false), pTparam0);
	//mum_PTadj = GausSmear(mum_PT, *pTparamsit);

	mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	musum = mup + mum;

	toys[0][toyindex]->Fill(musum.M());

	++toyindex;
      }

      if (npTmethods > 1) {
	toyindex=0;
	for ( auto pTparam1 : pTparams[1] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	  mup_PTadj = GausSmear_pdependent(GausSmear(mup_PT, mup_constsmear, false), pTparam1, mup_ETA);
	  //mup_PTadj = GausSmear_pTdependent(mup_PT, *pTparamsit);
	  mum_PTadj = GausSmear_pdependent(GausSmear(mum_PT, mum_constsmear, false), pTparam1, mum_ETA);
	  //mum_PTadj = GausSmear_pTdependent(mum_PT, *pTparamsit);

	  mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	  mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	  musum = mup + mum;

	  toys[1][toyindex]->Fill(musum.M());

	  ++toyindex;
	}
      }
      
      if (npTmethods > 2) {
	toyindex=0;
	for ( auto pTparam2 : pTparams[2] ) { //(pTparamsit = pTparams[2].begin(); pTparamsit != pTparams[2].end(); ++pTparamsit) {//
	  mup_PTadj = MomentumScale(mup_PTsmear, pTparam2);
	  //mup_PTadj = MomentumScale(mup_PTsmear, *pTparamsit);
	  mum_PTadj = MomentumScale(mum_PTsmear, pTparam2);
	  //mum_PTadj = MomentumScale(mum_PTsmear, *pTparamsit);

	  mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	  mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	  musum = mup + mum;

	  toys[2][toyindex]->Fill(musum.M());
	  
	  if (fill_asymsubsets) {
	    Double_t pasymmetry = (mum_PTadj*cosh(mum_ETA)-mup_PTadj*cosh(mup_ETA))/(mum_PTadj*cosh(mum_ETA) + mup_PTadj*cosh(mup_ETA));

	    if (pasymmetry < (-asym_limit)) {
	      asym_subsets[0][toyindex]->Fill(musum.M());
	    } else if  ((pasymmetry >= (-asym_limit)) && (pasymmetry < asym_limit)) {
	      asym_subsets[1][toyindex]->Fill(musum.M());
	    } else if (pasymmetry >= asym_limit) {
	      asym_subsets[2][toyindex]->Fill(musum.M());
	    } else {
	      cout << "Event " << eventit <<": Could not match the asymmetry to the desired range" << endl;
	    }
	  }
	  ++toyindex;
	}
      }
      
      if (npTmethods > 3) {
	toyindex=0;
	for ( auto pTparam3 : pTparams[3] ) { //(pTparamsit = pTparams[3].begin(); pTparamsit != pTparams[3].end(); ++pTparamsit) { //
	  //coindex = 0;
	  mup_PTadj = CurvatureBias(mup_PTsmear, pTparam3, echarge, mup_ETA);
	  mum_PTadj = CurvatureBias(mum_PTsmear, pTparam3, -echarge, mum_ETA);

	  mup.SetPtEtaPhiM(mup_PTadj, mup_ETA, mup_PHI, muMass);
	  mum.SetPtEtaPhiM(mum_PTadj, mum_ETA, mum_PHI, muMass);

	  musum = mup + mum;
	  toys[3][toyindex]->Fill(musum.M());
	  
	  // Fill histograms with desired observables for discerning curvature bias  
	  if (fill_curveobs || fill_profiles || fill_asymsubsets || fill_mscbtoys) {
	    Double_t pasymmetry = (mum_PTadj*cosh(mum_ETA)-mup_PTadj*cosh(mup_ETA))/(mum_PTadj*cosh(mum_ETA) + mup_PTadj*cosh(mup_ETA));
	    Double_t deltamup = (mum_PTadj*cosh(mum_ETA)-mup_PTadj*cosh(mup_ETA));
	    
	    if (fill_curveobs) {
	      curveobs[0][toyindex]->Fill(pasymmetry);
	      curveobs[1][toyindex]->Fill(deltamup);
	    }
	    if (fill_profiles) {
	      // We take the modulus in order to use only the odd numbered toys so that the graph will be less crowded
	      // ^ this functionality was moved into tcanvas_from_vector by introducing a step by which you skip some toys
	      //if ((toyindex % 2) == 1) {
		profiles[toyindex]->Fill(pasymmetry, musum.M());
		//}
	    }
	    if (fill_asymsubsets && (!fill_mscbtoys)) {
	      if (pasymmetry < (-asym_limit)) {
		asym_subsets[3][toyindex]->Fill(musum.M());
	      } else if  ((pasymmetry >= (-asym_limit)) && (pasymmetry < asym_limit)) {
		asym_subsets[4][toyindex]->Fill(musum.M());
	      } else if (pasymmetry >= asym_limit) {
		asym_subsets[5][toyindex]->Fill(musum.M());
	      } else {
		cout << "Event " << eventit <<": Could not match the asymmetry to the desired range" << endl;
	      }
	      // mscb = momentume scale curvature bias
	    } else if (fill_mscbtoys) {
	      mscbindex = 0;
	      for ( auto pscale : pTparams[2] ) {
		// mup_PTadj already has curvature bias applied, we simply add the momentum scaling
		mup.SetPtEtaPhiM(MomentumScale(mup_PTadj, pscale), mup_ETA, mup_PHI, muMass);
		mum.SetPtEtaPhiM(MomentumScale(mum_PTadj, pscale), mum_ETA, mum_PHI, muMass);
		musum = mup + mum;

		// Sory the results according to momentum asymmetry
		if (pasymmetry < (-asym_limit)) {
		  mscbtoys[mscbindex][toyindex][0]->Fill(musum.M());
		} else if  ((pasymmetry >= (-asym_limit)) && (pasymmetry < asym_limit)) {
		  mscbtoys[mscbindex][toyindex][1]->Fill(musum.M());
		} else if (pasymmetry >= asym_limit) {
		  mscbtoys[mscbindex][toyindex][2]->Fill(musum.M());
		} else {
		  cout << "Event " << eventit <<": Could not match the asymmetry to the desired range" << endl;
		}
		++mscbindex;
	      }
	    }
	  }
	  ++toyindex;
	}
      }
    } 

    if ((eventit % 2 == 1) || use_all_events) {
      //apply smear to template to better replicate experimental conditions - if smear_template=false smear value will just be the original
      mup.SetPtEtaPhiM(mup_PTsmear, mup_ETA, mup_PHI, muMass);
      mum.SetPtEtaPhiM(mum_PTsmear, mum_ETA, mum_PHI, muMass);

      musum = mup + mum;
      //cout << "Filling template histogram with invariant mass " << musum.M() <<endl;
      nominalH->Fill(musum.M());
    }
  }
  cout << "Histograms should all be filled: printing information" << endl;
  nominalH->Print("all");
  draw_colourhist2d(kinematics2d);

  for (vector< vector<TH1F*>>::iterator toyit = toys.begin(); toyit != toys.end(); ++toyit) {
    tcanvas_from_vector(*toyit,3);
  }
  if (fill_profiles) {
    tcanvas_from_vector(profiles,2,true,pTparams[3],false);
    //tcanvas_from_vector(profiles,2);
  }
  if (fill_asymsubsets) {
    tcanvas_from_vector(asym_subsets[0],3);
    tcanvas_from_vector(asym_subsets[1],3);
    tcanvas_from_vector(asym_subsets[2],3);
    tcanvas_from_vector(asym_subsets[3],3);
    tcanvas_from_vector(asym_subsets[4],3);
    tcanvas_from_vector(asym_subsets[5],3);
  } 

  ////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// PERFORMING CHI2 TEST /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  string chi2file = "./chi2results/" + propagator + "sampletest_lhcbcuts.csv";
  ofstream chi2_ofs (chi2file, ofstream::out);

  /*    if (data_it == pTparams[vect_row].begin()) {
	pT_ofs << *data_it; 
      } else {
	pT_ofs << ',' << *data_it; 
      }
    }
    pT_ofs << endl;  
  }
  
  pT_ofs.close();
  */
  
  Double_t chi2point=0;
  // variable rechi2 will be the test statistic for the toy with all neutral parameters, which we will use when we make a 2d histogram
  // of the chi2 for different momentum scale and curvature bias values
  Double_t refchi2;
  //Double_t event_count_exp = lhcb_luminosity*xs_Zmumu;
  Double_t template_count, toy_count;
  
  for (int pTmethod = 0; pTmethod < npTmethods; ++pTmethod) {
    toyindex=0;
    for (toyit = toys[pTmethod].begin(); toyit != toys[pTmethod].end(); toyit++) {
      if((*toyit)->GetSumw2N() == 0) {
	cout << "Warning! Weights do not seem to be stored" << endl;
      }

     template_count = nominalH->Integral();
     toy_count = (*toyit)->Integral();
     //cout << "Scaling template integral " << template_count << " and toy " << toy_count << " to expected count" << endl;
     nominalH->Scale(toy_count / template_count);       
     nbins = nominalH->GetNbinsX();
       
     if (use_all_events) {   
       for (int binit = 0; binit < nbins; ++binit) {
	 //Double_t template_bin = nominalH->GetBinContent(binit);
	 //Double_t template_error = nominalH->GetBinError(binit);
	 Double_t toy_bin = (*toyit)->GetBinContent(binit);
      
	 (*toyit)->SetBinContent(binit, gRandom->Poisson(toy_bin));
	 toy_bin = (*toyit)->GetBinContent(binit);

	 Double_t toy_error = sqrt(toy_bin);
	 (*toyit)->SetBinError(binit, toy_error);
	 //chi2point += (pow(toy_bin-template_bin,2))/(pow(toy_error,2) + pow(template_error, 2));
       }
       chi2point = (*toyit)->Chi2Test(nominalH, "Chi2 WW");
     
     } else {
       chi2point = (*toyit)->Chi2Test(nominalH, "Chi2 WW");
       // set refchi2 at index 5 because this is the central value (1 for momentum scale, 0 for curvature bias)
       if (toyindex==5) {
	 refchi2 = chi2point;
	 cout << "Reference chi2 test statistic is: " << refchi2 << endl;
       }
     }
     
     cout << "Adding method " << pTmethod << " toy " << toyindex << " result to csv: " << chi2point << endl;
     if (toyindex == 0) {
       chi2_ofs << chi2point;
     } else {
       chi2_ofs << ',' << chi2point;
     }
     ++toyindex;
    }
    chi2_ofs << endl;
  }

  if (fill_curveobs){
    for (coindex = 0; coindex < 2; ++coindex) {
      toyindex=0;
      for (toyit = next(curveobs[coindex].begin()); toyit != curveobs[coindex].end(); ++toyit) {
	chi2point = (*toyit)->Chi2Test((*curveobs[coindex].begin()), "Chi2 WW");

	if (toyindex == 0) {
	  chi2_ofs << chi2point;
	} else {
	  chi2_ofs << ',' << chi2point;
	}
	++toyindex;
      }
      chi2_ofs << endl;
    }
  }

  Double_t chi2_asym;
  if (fill_mscbtoys) {
    //Double_t stepms = (pTparams[2][-1]-pTparams[2][0])/(pTparams[2].size() -1);
    //Double_t stepcb = (pTparams[3][-1]-pTparams[3][0])/(pTparams[3].size() -1);
    //TH2F *mscb_chi2hist = new TH2F("mscb_chi2hist", "chi2 comparison with nominal template", 6, pTparams[2][0], pTparams[2][-1] + stepms, 6, pTparams[3][0], pTparams[3][-1] + stepcb);
    Double_t xedges[12] = {0.999934, 0.999946, 0.999958, 0.999970, 0.999982, 0.999994, 1.000006, 1.000018, 1.000030, 1.000042, 1.000054, 1.000066};
    Double_t yedges[12] = {-1.1e-25, -0.9e-25, -0.7e-25, -0.5e-25, -0.3e-25, -0.1e-25, 0.1e-25, 0.3e-25, 0.5e-25, 0.7e-25, 0.9e-25, 1.1e-25};
    char mscbhist_name[30];
    sprintf(mscbhist_name, "mscb_chi2hist%s", propagator.c_str());
    TH2F *mscb_chi2hist = new TH2F(mscbhist_name, "chi2 comparison with nominal template", 11, xedges, 11, yedges);
    mscb_chi2hist->GetXaxis()->SetTitle("Momentum scale");
    mscb_chi2hist->GetYaxis()->SetTitle("Curvature Bias");

    for (int msindex = 0; msindex < pTparams[2].size(); ++msindex) {
      for (int cbindex = 0; cbindex < pTparams[3].size(); ++cbindex) {
	chi2point = 0;
												    
	for (int asindex = 0; asindex < 3; ++asindex) {
	  TH1F *mscbtoy = new TH1F(*mscbtoys[msindex][cbindex][asindex]);
	  nominalH->Scale( mscbtoy->Integral() / nominalH->Integral() );
	  chi2_asym = mscbtoy->Chi2Test(nominalH, "Chi2 WW");
	  cout << "MS = " << pTparams[2][msindex] << ", CB = " << pTparams[3][cbindex] << ", asym subset = " << asindex<< ", CHI2 = " <<  chi2_asym << endl;
	  
	  //We take the sum across 3 asym subsets for the chi2 test because curvature bias not detectable in combined results
	  // or in the samples with asymmetry near zero.  This way elevated chi-square in each of the extreme subsets will have 
	  // an effect on the total	
	  chi2point += chi2_asym;
	}
	// Take the logarithm of the difference between chi2 with these parameters and the reference
	mscb_chi2hist->SetBinContent(msindex+1, cbindex+1, log10(chi2point - refchi2));
      }
    }
    draw_colourhist2d(mscb_chi2hist);
  }


  chi2_ofs.close();
  output->Write();
  output->Close();
}

/*
void Zanalysis (string pTparamfile="./pTparameters.csv", string rootfile="/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root") {

  vector < vector<Double_t>> pTparams;
  ////////////// PUT CSV READING CODE HERE /////////////////////
  
  ifstream input;
  input.open(pTparamfile);

  if (input.fail()) {
    cerr << "Error in pT parameter file" <<endl;
    exit(1)
  } 
   
  ///////// from here assume pT params is initialized for the time being ////////////////////////////////////////////
  if (pTparams.size() < 1) {
    return;
  } else {
    Zanalysis(pTparams, rootfile);
  }
}
*/


void Wsample_analysis (int ndata_files, string Wcharge,string dimu_propagator,
		       string dimurootfile="",
		       bool smear_template = true,
		       bool saveWroot = true,
		       bool fill_curveobs = false,
		       bool fill_profiles = false,
		       bool fill_asymsubsets = false,
		       bool fill_mscbtoys = false,
		       bool fill_smeartests = false,
		       bool dimuuse_all_events = false) {

  int hist_dims[3] = {40,30,50};

  //int ndata_files = 20;
  //string Wcharge = "Wm";
  string pTmethods[4] = {"GausSmear", "GausSmear_pTdependent", "MomentumScale", "CurvatureBias"};
  int npTmethods = 4;

  int ntemplates = 11;
  int ntoys = 11;
  double pTparam_limits[8] = {0.0, 0.004, 0.0, 0.0001, 0.99994, 1.00006, -1e-25, 1e-25};

  Double_t MWnom = 80.4;  
  Double_t gamma = 2.15553;
  Double_t echarge = 1.602176565e-19;

  Double_t lhcb_luminosity = 6.0; //Predicted fb^-1 at the end of run 2
  Double_t xs_Wp = 1093600.0, xs_Wm = 818400; //quantities in fb
  bool use_all_events = false;

  TH1::SetDefaultSumw2();
  vector<TH1F *> templates;
  vector< vector<TH1F *> > toys(npTmethods);
  vector<Double_t> Wmasses;
  vector< vector<Double_t> > pTparams(npTmethods);

  vector<TH1F *>::iterator toyit, templateit;
  vector<Double_t>::iterator Wmassit, pTparamsit;

  string tbuffer = make_timestamp();

  string fileinfo = Wcharge + "_" + to_string(ndata_files) + "ntuples_lhcbcuts_" + tbuffer;
  string output_name = "~/oxford-lhcb-wmass/rootfiles/" + fileinfo + ".root";

  //create output file  	
  TFile *output = new TFile(output_name.c_str(),"RECREATE");

  string Wnominalstr;
  Wnominalstr = trunc_double(MWnom, 3); 
  int toyindex=0, templateindex=0;
  string template_name;
  string template_title;

  // Initialize a vector of W mass hypotheses which will assist in iterating over and filling template histograms
  for (Double_t MWhyp=79.8; MWhyp<=80.8; MWhyp+=(80.8-79.8)/(ntemplates-1)) { //Mhyp for mass hypothesis
    Wmasses.push_back(MWhyp);  

    //Wreweightss << fixed <<setprecision(3) << MWhyp;
    template_name = Wcharge + "template" + to_string(templateindex);
    template_title = "mu_PT with W mass reweight " + trunc_double(MWhyp, 3) + " - nominal " + Wnominalstr; 

    TH1F *hweighted_template = new TH1F(template_name.c_str(), template_title.c_str(), hist_dims[0], hist_dims[1], hist_dims[2]);    
    templates.push_back(hweighted_template);
   
    //Resetting the stringstream
    //Wreweightss.str(string());
    ++templateindex;
  }
  templateindex = 0;
  
  //stringstream value;
  string toy_name;
  char toy_title[100];

  Double_t step0, step1, step2, step3;
  step0 = (pTparam_limits[1]-pTparam_limits[0])/(ntoys-1);
  step1 = (pTparam_limits[3]-pTparam_limits[2])/(ntoys-1);
  step2 = (pTparam_limits[5]-pTparam_limits[4])/(ntoys-1);
  step3 = (pTparam_limits[7]-pTparam_limits[6])/(ntoys-1);

  for (Double_t pTparam0 = pTparam_limits[0]; pTparam0 <= (pTparam_limits[1] + 0.5*step0); pTparam0 += step0) {
    pTparams[0].push_back(pTparam0);
    //value << fixed << setprecision(3) << pTparam;
    toy_name =  pTmethods[0]  + Wcharge + to_string(toyindex);
    sprintf(toy_title, "sigma mu_PT ~ %f", pTparam0);
    TH1F *hpT_toy = new TH1F(toy_name.c_str(), toy_title, hist_dims[0], hist_dims[1], hist_dims[2]);
    toys[0].push_back(hpT_toy);
    
    ++toyindex;
  }

  if (npTmethods > 1) {
    toyindex = 0;
    for (Double_t pTparam1 = pTparam_limits[2]; pTparam1 <= (pTparam_limits[3] + 0.5*step1); pTparam1 += step1) {
      pTparams[1].push_back(pTparam1); 

      //value << fixed << setprecision(5) << pTparam;
      toy_name = pTmethods[1]  + Wcharge + to_string(toyindex);
      sprintf(toy_title, "sigma mu_PT ~ %f * mu_PT", pTparam1);
      TH1F *hpT_toy = new TH1F(toy_name.c_str(), toy_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[1].push_back(hpT_toy);

      ++toyindex;
    }
  }

  if (npTmethods > 2) {
    toyindex = 0;
    for (Double_t pTparam2 = pTparam_limits[4]; pTparam2 <= (pTparam_limits[5] + 0.5*step2); pTparam2 += step2) {
      pTparams[2].push_back(pTparam2);
            //value << fixed << setprecision(4) << pTparam;
      toy_name = pTmethods[2] + Wcharge + to_string(toyindex);
      sprintf(toy_title, "mu_PT -> %f * mu_PT", pTparam2);
      TH1F *hpT_toy = new TH1F(toy_name.c_str(), toy_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[2].push_back(hpT_toy);

      ++toyindex;
    }
  }

  if (npTmethods > 3) {
    toyindex=0;
    for (Double_t pTparam3 = pTparam_limits[6]; pTparam3 <= (pTparam_limits[7] + 0.5*step3); pTparam3 += step3) {
      pTparams[3].push_back(pTparam3);

      toy_name = pTmethods[3] + Wcharge + to_string(toyindex);
      sprintf(toy_title, "charge/mu_PT -> charge/mu_PT + %s", trunc_double(pTparam3, 2, "scientific").c_str());
      TH1F *hpT_toy = new TH1F(toy_name.c_str(), toy_title, hist_dims[0], hist_dims[1], hist_dims[2]);
      toys[3].push_back(hpT_toy);

      ++toyindex;
    }
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
  Float_t prop_M, mu_PT, mu_ETA;
  Float_t mu_PTadj, mu_PTsmear, mu_PTsmeartest;
  MCDecayTree->SetBranchAddress("mu_PT", &mu_PT);
  MCDecayTree->SetBranchAddress("mu_ETA", &mu_ETA);
  MCDecayTree->SetBranchAddress("prop_M", &prop_M);

  MCDecayTree->SetBranchStatus("*", 0);
  MCDecayTree->SetBranchStatus("mu_PT", 1);
  MCDecayTree->SetBranchStatus("mu_ETA", 1);
  MCDecayTree->SetBranchStatus("prop_M",1);

  Long64_t nentries = MCDecayTree->GetEntries();
  cout << "nentries: " << nentries << endl;
  Long64_t nbytes = 0;
  Double_t constsmear, psmear;

  ////////////////// Initialize histograms of the total momentum and a Profile of resolution  versus momentum /////////
  TH1F *momentumH = new TH1F("momentumH", "Total momentum of muon from W decay", 50, 50, 500);
  //TProfile *dppProfile = new TProfile("dppProfile", "Profile of dp/p versus p", 50, 0, 500);

  /////////// Initialize any 2D histograms /////////////////
  string kinematics2d_name = Wcharge + "kin2d_logpT_eta";
  // Note that specific limits need to be specified for both axes for it to work (i.e. you cannot do automatic axis range for one and not for the other)
  TH2F *kinematics2d = new TH2F(kinematics2d_name.c_str(), "Histogram of pT and eta of muon from W decay", 15, 0, 1, 10, 2.0, 4.5);
  kinematics2d->GetXaxis()->SetTitle("log10(pT) (GeV/c)");
  kinematics2d->GetYaxis()->SetTitle("eta");

  /////////////// INITIALIZE HISTOGRAMS FOR VERIFYING EQUIVALENCE OF ADDING WEIGHTED SMEAR COEFFICIENTS AND DRAWING FROM GAUSSIAN WITH PROPAGATED UNCERTAINTY VARIANCE ////////////////
  
  //if (fill_smeartests){
    vector<TH1F *> smeartests(2);
    smeartests[0] = new TH1F("smeartest1", "Constant smear then pT dependent smear", 40, 30, 50);
    smeartests[1] = new TH1F("smeartest2", "pT dependent smear then constant smear", 40, 30, 50);
    //}

  //Get the correct cross section for luminosity scaling of chi square (fcurrently not implemented)
  Double_t xs;
  if(Wcharge.compare("Wm") == 0) {   
    echarge = -echarge;
    xs = xs_Wm;
  } else {
    xs = xs_Wp;
  }

  for (Long64_t eventit=0; eventit < nentries; eventit++) { //usually i< nentries*split_ratio for full data set
    nbytes += MCDecayTree->GetEntry(eventit);

    if (smear_template) {
      constsmear = gRandom->Gaus(1, 0.5*(pTparams[0].back()));      
      psmear = gRandom->Gaus(1, 0.5*(pTparams[1].back())*mu_PT*cosh(mu_ETA));      
      mu_PTsmear = GausSmear(mu_PT, constsmear, false);
      mu_PTsmear = GausSmear_pdependent(mu_PTsmear, psmear, mu_ETA, false);
      
      if (fill_smeartests) {
	smeartests[0]->Fill(mu_PTsmear);
	//The intention here was to test if changing the order of these smears would impact the shape of the histograms
	smeartests[1]->Fill(GausSmear(GausSmear_pdependent(mu_PT, 0.5*pTparam_limits[3], mu_ETA), 0.5*pTparam_limits[1], mu_ETA));      
      }     
    } else {
      mu_PTsmear = mu_PT;
    }
      
    // Find events which do not meet the LHCb acceptance range or minimu pT and skip them
    if ((mu_PTsmear < 20) || (mu_ETA < 2.0) || (mu_ETA > 4.5)) {
      continue;
    }
    momentumH->Fill(mu_PT*cosh(mu_ETA));
    kinematics2d->Fill(log10(mu_PT), mu_ETA);
    //propM->Fill(prop_M);

    if ((eventit%2 == 0) || use_all_events) {
      
      toyindex=0;
      for ( auto pTparam0 : pTparams[0] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	mu_PTadj = GausSmear_pdependent(mu_PT, psmear, mu_ETA, false);
	mu_PTadj = GausSmear(mu_PTadj, pTparam0, mu_PT);
	
	toys[0][toyindex]->Fill(mu_PTadj);
	++toyindex;
      }
      
      if (npTmethods > 1) {
	toyindex=0;
	for ( auto pTparam1 : pTparams[1] ) {  // (pTparamsit = pTparams.begin(); pTparamsit != pTparams.end(); ++pTparamsit) {
	  mu_PTadj = GausSmear(mu_PT, constsmear, false);
	  mu_PTadj = GausSmear_pdependent(mu_PTadj, pTparam1, mu_ETA);
	  //cout << "mu pT after curvature resolution adjustment " << mu_PTadj << endl;
	  toys[1][toyindex]->Fill(mu_PTadj);
	  ++toyindex;
	}
      }
	
      if (npTmethods > 2) {
	toyindex=0;
	for ( auto pTparam2 : pTparams[2] ) {
	  mu_PTadj = MomentumScale(mu_PTsmear, pTparam2);
	  toys[2][toyindex]->Fill(mu_PTadj);
	  ++toyindex;
	}
      }
	
      if (npTmethods > 3) {
	toyindex=0;
	for ( auto pTparam3 : pTparams[3] ) {
	  mu_PTadj = CurvatureBias(mu_PTsmear, pTparam3, echarge, mu_ETA);
	  toys[3][toyindex]->Fill(mu_PTadj);
	  ++toyindex;
	}
      }
    }  

    if ((eventit % 2 == 1) || use_all_events) {
      templateindex = 0;
      const auto denominator = TMath::BreitWigner(prop_M, MWnom, gamma); 
      //mu_PTadj = GausSmear(mu_PTsmear, 0.005);
      //mu_PTadj = mu_PT;
      for ( auto Wmass : Wmasses ) {
	//cout << "Filling template histogram for W mass " << Wmass << endl;
	templates[templateindex]->Fill(mu_PTsmear, TMath::BreitWigner(prop_M, Wmass, gamma)/denominator);
	++templateindex;
      }
    }
  }
 
  string templatesHname = "~/oxford-lhcb-wmass/plots/" + Wcharge + "templatesHist.pdf";
  tcanvas_from_vector(templates, 2, true, Wmasses,"cBW", "mu PT with different W mass hypotheses", templatesHname);
  for (vector< vector<TH1F*>>::iterator toyit = toys.begin(); toyit != toys.end(); ++toyit) {
    tcanvas_from_vector(*toyit,2);
  }
	 
  if (fill_smeartests) {
    smeartests[0]->Scale(1/(smeartests[0]->Integral()));
    smeartests[1]->Scale(1/(smeartests[1]->Integral()));
    tcanvas_from_vector(smeartests);
  }
  draw_colourhist2d(kinematics2d);
  /*
  TCanvas *kin2dc = new TCanvas("kin2dc", kinematics2d->GetTitle());
  gStyle->SetPalette(57);
  kinematics2d->SetContour(20);
  kinematics2d->Draw("colz");
  kin2dc->Print("/home/egan/oxford-lhcb-wmass/Wkinematics2d.pdf");
  kin2dc->Print("/home/egan/oxford-lhcb-wmass/Wkinematics2d.png");
  kin2dc->Close();
  */

  /////////////////////////////// PERFORMING A CHI SQUARE TEST //////////////////////////////////////////

  cout << "All vectors and iterators initialized.  template_vect is of length " << templates.size() << ". toy_vect is of length " <<  toys.size() << ". Wmass_vect is of length " << Wmasses.size() << endl;
  
  TGraph *chi2Plot = new TGraph(ntemplates);

  char chi2plot_name[50];
  Double_t chi2point=0;
  toyindex=0;
  templateindex=0;  

  Double_t template_count, toy_count;
  int nbins; 

  TCanvas* c2 = new TCanvas("cchi2", "chi2 with different W mass hypotheses");
  
  for (int pTmethod = 0; pTmethod < npTmethods; ++pTmethod) {
    toyindex=0;
    for (toyit = toys[pTmethod].begin(); toyit != toys[pTmethod].end(); toyit++) {

      if ((*toyit)->GetSumw2N() == 0) {
       cout << "Warning! Weights do not seem to be stored" << endl;
     }

     for (templateit = templates.begin(); templateit != templates.end(); templateit++){
       
       if (((*toyit)->Integral() != 0) && ((*templateit)->Integral() != 0)) {
	 (*templateit)->Scale((*toyit)->Integral()/(*templateit)->Integral());
       } else {
	 break;
       }
              
       nbins = (*templateit)->GetNbinsX();       
       sprintf(chi2plot_name, "chi2plot%s%u%u", Wcharge.c_str(), pTmethod, toyindex);
       
       // METHOD 1: Go bin by bin summing the squared difference of bins divided by the sum of errors squared
       if (use_all_events) {
	 
	 for (int binit = 0; binit < nbins; ++binit) {
	   Double_t template_bin = (*templateit)->GetBinContent(binit);
	   Double_t template_error = (*templateit)->GetBinError(binit);
	   Double_t toy_bin = (*toyit)->GetBinContent(binit);
	 
	   (*toyit)->SetBinContent(binit, gRandom->Poisson(toy_bin));
	   toy_bin = (*toyit)->GetBinContent(binit);

	   Double_t toy_error = sqrt(toy_bin);
	   (*toyit)->SetBinError(binit, toy_error);
	   Double_t chi2num = pow(toy_bin - template_bin, 2);
	   Double_t chi2denom = pow(toy_error, 2) + pow(template_error, 2);
	   chi2point += chi2num/chi2denom;
	   cout << "Bin " << binit << ": adding to chi2 sum  " << chi2num << "/" << chi2denom << endl;
	   cout << "Unreduced chi2 statistic: " << chi2point << endl;
	 }
	 chi2point = chi2point/(nbins-2);

       } else {
	 chi2point = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
       }
       // METHOD 2: apply root chi2 test function - returning chi2 test statistic (as opposed to p-value) for weighted comparison
       //chi2point1 = (*toyit)->Chi2Test((*templateit), "Chi2 WW");
       
       cout << "Copying chi2 result 1 for toy " << toyindex << " and template " << templateindex << ": " << chi2point << endl;

       chi2Plot->SetPoint(templateindex, Wmasses[templateindex], chi2point);

       templateindex = (templateindex + 1) % templates.size();
     }

     cout << "Drawing graph of toy " << toyindex << endl;
     chi2Plot->SetMarkerColor(4);
     
     chi2Plot->SetMarkerStyle(8);
     
     chi2Plot->Draw("AP");
     
     cout << "Writing chi2 plot to root file." << endl;
     output->WriteTObject(chi2Plot, chi2plot_name, "Overwrite");
     cout << "Write successful" << endl;

     //c2->Print("./plots/chi2Plot.png");
     c2->Clear();  

     toyindex++;
   }
  }
  c2->Close();
  
  if (saveWroot) {	 
    output->Write();
  }
  output->Close();

  if ((dimu_propagator.compare("Y") == 0) || (dimu_propagator.compare("Z") == 0)) {
    dimuon_analysis(pTparams, dimu_propagator, dimurootfile, smear_template, fill_curveobs, fill_profiles, fill_asymsubsets, fill_mscbtoys, dimuuse_all_events);
  }
}

int main (int argc, const char** argv) {
  int ndata_files;
  string Wcharge, Zrootfile, Yrootfile, dimu_propagator;
  bool saveWroot, dimuuse_all_events, smear_template, fill_curveobs, fill_profiles, fill_asymsubsets, fill_mscbtoys, fill_smeartests;

  po::options_description desc("Allowed options");
  
  desc.add_options()
    ("help", "produce help message")
    //("rootfile", po::value<string>(&rootfile)->default_value("/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root"))
    ("ndata_files", po::value(&ndata_files)->default_value(20))
    //("pTparamfile", po::value<string>(&pTparamfile)->default_value("./pTparameters/pTparameters.csv"))
    ("Wcharge", po::value(&Wcharge)->default_value("Wp"))
    ("saveWroot", po::value(&saveWroot)->default_value(true))
    //("Zoutput_name", po::value(&Zoutput_name)->default_value("./rootfiles/Zsampletest_lhcbcuts.root"))
    ("Zrootfile", po::value(&Zrootfile)->default_value(""))
    ("Yrootfile", po::value(&Yrootfile)->default_value(""))
    ("dimu_propagator", po::value(&dimu_propagator))
    ("dimuuse_all_events", po::value(&dimuuse_all_events)->default_value(false))
    ("smear_template", po::value(&smear_template)->default_value(true))
    ("fill_curveobs", po::value(&fill_curveobs)->default_value(false))
    ("fill_asymsubsets", po::value(&fill_asymsubsets)->default_value(false))
    ("fill_mscbtoys", po::value(&fill_mscbtoys)->default_value(false))
    ("fill_profiles", po::value(&fill_profiles)->default_value(false))
    ("fill_smeartests", po::value(&fill_smeartests)->default_value(false))
  ;
  /*
  po::positional_options_description pod;
  pod.add("ndata_files", -2);
  pod.add("Wcharge", -1);
  */

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);


  if (fill_mscbtoys && fill_asymsubsets) {
    cout << "WARNING: simultaneous filling of mscbtoys and asymsubsets is currently not supported.  Note that mscbtoys contains all the info in asymsubsets, but in a way that is more difficult to extract" << endl;
  }
  if (!(vm.count("dimu_propagator"))) {
    cout << "WARNING: must specify dimu propagator" << endl;
  } else if (dimu_propagator.compare("Z") == 0) {
    Wsample_analysis(ndata_files, Wcharge, dimu_propagator, Zrootfile, smear_template, saveWroot,
		     fill_curveobs, fill_profiles, fill_asymsubsets, fill_mscbtoys, fill_smeartests, dimuuse_all_events);
  } else if (dimu_propagator.compare("Y") == 0) {
    Wsample_analysis(ndata_files, Wcharge, dimu_propagator, Yrootfile, smear_template, saveWroot,
		     fill_curveobs, fill_profiles, fill_asymsubsets, fill_mscbtoys, fill_smeartests, dimuuse_all_events);
  } else {
    Wsample_analysis(ndata_files, Wcharge, dimu_propagator, "", smear_template);
  }
  return 0;
}
