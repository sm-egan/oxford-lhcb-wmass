#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLorentzRotation.h"
#include "TSystem.h"
#include <vector>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>
#include <typeinfo>

///data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root

namespace po = boost::program_options;
using namespace std;


int main ( int argc, char* argv[]) {

  string rootfile, pTparamfile;
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    //("rootfile", po::value<string>(&rootfile)->default_value("/data/lhcb/users/mvesteri/GenLevelV19/merged/13TeV_Z_PowhegPythiaDefault.root"))
    ("rootfile", po::value< vector<string> >())//&rootfile)->default_value("./rootfiles/Zsampletest.root"))
    //("pTparamfile", po::value<string>(&pTparamfile)->default_value("./pTparameters/pTparameters.csv"))
    ("pTparamfile", po::value< vector<string> >())//&pTparamfile)->default_value("./pTparameters/pTparameters.csv"))
  ;

  po::positional_options_description p;
  p.add("rootfile", -2);
  p.add("pTparamfile", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  
  cout << "CHECKING THAT COMMAND LINE OPTIONS HAVE BEEN PARSED PROPERLY" << endl;
  
  /*
  if (vm.count("rootfile")) {
    cout << "Rootfile was set to " << vm["rootfile"].as< vector<string> > << " and has type" << typeid(vm["rootfile"]).name() << endl;
    
  } else {
    cout << "rootfile was not set";
  }

  if (vm.count("pTparamfile")) {
    cout << "Rootfile was set to " << vm["pTparamfile"].as< vector<string> > << " and has type" << typeid(vm["pTparamfile"]).name() << endl;
    
  } else {
    cout << "pT parameter file was not set";
  }
  */

  cout << "rootfile is set to: " << rootfile << endl;
  cout << "pTparamfile is set to: " << pTparamfile << endl;
  
//Zanalysis(pTparamfile, rootfile);
}
