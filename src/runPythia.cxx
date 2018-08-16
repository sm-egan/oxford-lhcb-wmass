#include "Pythia8/Pythia.h"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main( int argc, const char** argv )
{
  
  int nEvents, seed;
  std::string process,MPI;
  std::string outputFile;
  float mHatMin,mHatMax,eCM;
  std::string alphaS,IKT;//to avoid truncation mismatch between file name and true used value!
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("nEvents",po::value(&nEvents)->default_value(1000))
    ("seed",po::value(&seed)->default_value(1))
    ("eCM",po::value(&eCM)->default_value(13000.)) // COM enery 13 TeV
    ("outputFile",po::value(&outputFile)->default_value("/home/egan/oxford-lhcb-wmass/Ysamples/pythia_upsilon.root"))
    ;

  
  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  po::notify( vm );

  if ( vm.count( "help" ) ) {
    std::cout << desc << std::endl;
    return 1;
  }
  
  Pythia8::Pythia pythia;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Main:numberOfEvents = "+std::to_string(nEvents));
  pythia.readString("Beams:eCM = "+std::to_string(eCM));
  pythia.readString("Random:seed = "+std::to_string(seed));
  pythia.readString("Bottomonium:all = off");//TODO: figure out how to switch on only Upsilon(1S)
  pythia.readString("Bottomonium:gg2bbbar(3S1)[3S1(1)]g = on,off,off");
  pythia.readString("553:onMode = off");
  pythia.readString("553:onIfAny = 13"); //only decays to muon for the Upsilon(1S) (MC particle ID = 553)
  pythia.readString("PartonLevel:MPI = off");

  TFile * output = new TFile(outputFile.c_str(),"RECREATE");
  TTree * tree = new TTree("MCDecayTree","MCDecayTree");

  Float_t prop_M,
    lep1_PT,
    lep1_ETA,
    lep1_PHI,
    lep1_ID,
    lep2_PT,
    lep2_ETA,
    lep2_PHI,
    lep2_ID,
    dilep_M,
    dilep_PT,
    dilep_Y;
  

  const TString lep1brname ("mup");
  tree->Branch(lep1brname +"_PT",&lep1_PT,lep1brname +"_PT/F");
  tree->Branch(lep1brname +"_ETA",&lep1_ETA,lep1brname +"_ETA/F");
  tree->Branch(lep1brname +"_PHI",&lep1_PHI,lep1brname +"_PHI/F");
  tree->Branch(lep1brname +"_ID",&lep1_ID,lep1brname +"_ID/I");

  const TString lep2brname ("mum");
  tree->Branch(lep2brname +"_PT",&lep2_PT,lep2brname +"_PT/F");
  tree->Branch(lep2brname +"_ETA",&lep2_ETA,lep2brname +"_ETA/F");
  tree->Branch(lep2brname +"_PHI",&lep2_PHI,lep2brname +"_PHI/F");
  tree->Branch(lep2brname +"_ID",&lep2_ID,lep2brname +"_ID/I");
    
  pythia.init();
  
  TLorentzVector lep1,lep2;
  std::map <double, Pythia8::Particle> mu_map;
  std::map <double, Pythia8::Particle> antimu_map;
  std::pair <std::map<double, Pythia8::Particle>::iterator,bool> inserted;

  for (int iEvent = 0; iEvent < nEvents; ++iEvent ) {
    antimu_map.clear();
    mu_map.clear();
    //std::cout << "////////////////////// event number " << iEvent << " //////////////////////" << std::endl;
    if (!pythia.next()) continue;    
    
    for ( auto i = 0; i < (int)pythia.event.size(); ++i){
      Pythia8::Particle P(pythia.event[i]);
      if(P.isFinal()){
	bool goodLep1 = P.id() == -13; // look specifically for antimuon
	bool goodLep2 = P.id() == 13; // look specifically for muon
	if (goodLep1) {//&& (P.pT() > lep1.Pt())) {
	  
	  inserted = antimu_map.insert( std::pair<double, Pythia8::Particle>(P.pT(), P) );
	  /*
	  if (inserted.second == false) {
	    std::cout << "Antimuon with identicle pT " << P.pT() << " was previously inserted in the map - keys are degenerate" << std::endl;
	  } else {
	    std::cout << "ADDING TO MAP ANTIMU WITH PT: " << P.pT() << std::endl;
	    }*/
	  
	} else if (goodLep2){ //&& (P.pT() > lep2.Pt())){

	  inserted = mu_map.insert( std::pair<double, Pythia8::Particle>(P.pT(), P) );
	  /*
	  if (inserted.second == false) {
	    std::cout << "Muon with identicle pT " << P.pT() << " was previously inserted in the map - keys are degenerate" << std::endl;
	  } else {
	    std::cout << "ADDING TO MAP MU WITH PT: " << P.pT() << std::endl;
	    } */
	} 
      }
    }

    if (antimu_map.empty() || mu_map.empty()) {
      //std::cout << "Either muon or antimuon map is empty! Skip to next event" << std::endl;
      continue;
    } else {
      Pythia8::Particle antimu = Pythia8::Particle((*antimu_map.rbegin()).second);
      Pythia8::Particle mu = Pythia8::Particle((*mu_map.rbegin()).second);
      //std::cout << "ADDING TO TREE ANTIMU WITH PT: " << antimu.pT() << std::endl;
      //std::cout << "ADDING TO TREE MU WITH PT: " << mu.pT() << std::endl;
      lep1.SetPxPyPzE(antimu.px(), antimu.py(), antimu.pz(), antimu.e());
      lep2.SetPxPyPzE(mu.px(), mu.py(), mu.pz(), mu.e());
      
      lep1_PT = lep1.Pt();
      lep1_ETA = lep1.Eta();
      lep1_PHI = lep1.Phi();
      lep2_PT = lep2.Pt();
      lep2_ETA = lep2.Eta();
      lep2_PHI = lep2.Phi();
    
      TLorentzVector lepsum = lep1+lep2;
      //std::cout << "INVARIANT MASS OF MUON SYSTEM IS: " << lepsum.M() << std::endl;
      
      tree->Fill();
    }
    
  } // End of event loop.
  
  output->Write();
  output->Close();
  pythia.stat();

  return 0;
}
