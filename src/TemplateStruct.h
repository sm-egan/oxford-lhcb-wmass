#ifndef TemplateStruct_h
#define TemplateStruct_h

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



using namespace std;

class TemplateStruct {

	public: 
		vector<Double_t> Wmasses; 
		vector<TH1F *> templates; 
		vector<TH1F *> toys;
		Double_t split_ratio;
		TString output_name;

	        TemplateStruct();
		TemplateStruct(Double_t split, TString root_file);
		
		void create_templates(int hist_dims[3], int ntemplates=5, int ndata_files=2);

	protected:
		void drawH_BW (TChain* EventChain, string ReweightBranch, string HistBranch, 
			       int hist_dims[3], Double_t nominal_mean, Double_t reweight_mean);

};

#endif
