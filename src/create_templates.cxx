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
  TemplateStruct *ts = new TemplateStruct();

  ts->create_templates(hist_dims, 11, 2);
  cout << "template_vect is of length " << ts->templates.size() << ". toy_vect is of length " <<  ts->toys.size() << ". Wmass_vect is of length " << ts->Wmasses.size() << endl;
  cout << (ts->templates)[0] << endl;
  cout << (ts->toys).front() << endl;
  //cout << "templates created. Now performing chi2 test" << endl;
  //ts->template_chi2();
}
