#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"

#include "ShipMCTrack.h"

// Cuts
#include "sndBaseCut.h"
#include "sndMinSciFiHitsCut.h"
#include "sndSciFiStationCut.h"
#include "sndVetoCut.h"
#include "sndMinSciFiConsecutivePlanes.h"
#include "sndDSActivityCut.h"
#include "sndUSQDCCut.h"
#include "sndEventDeltat.h"
#include "sndAvgSciFiFiducialCut.h"
#include "sndAvgDSFiducialCut.h"

int main(int argc, char ** argv) {

  std::cout << "Starting neutrino filter" << std::endl;
  
  if (argc != 3) {
    std::cout << "Two arguments required: input file name (or reg exp), output file name" << std::endl;
    exit(-1);
  }

  // Input files
  bool isMC = false;
  TChain * ch = new TChain("rawConv");
  ch->Add(argv[1]);
  if (ch->GetEntries() == 0){
    delete ch;
    ch = new TChain("cbmsim");
    if (ch->GetEntries() > 0) isMC = true;
    else {
      std::cout << "Didn't find rawConv or cbmsim in input file" << std::endl;
      exit(-1);
    }
  }
  std::cout << "Got input tree" << std::endl;

  // MC truth
  TClonesArray * MCTracks = new TClonesArray("ShipMCTrack", 5000);
  if (isMC) ch->SetBranchAddress("MCTrack", &MCTracks);

  // Output file
  TFile * outFile = new TFile(argv[2], "RECREATE");
  std::cout << "Got output file" << std::endl;  

  // Set up all branches to copy to output TTree.
  TTree * outTree = ch->CloneTree(0);
  std::cout << "Got output tree" << std::endl;  

  // Set up cuts
  std::cout << "Starting cut set up" << std::endl;
  std::vector< sndAnalysis::baseCut * > cutFlow;
  cutFlow.push_back( new sndAnalysis::minSciFiHits(1, ch)); // A. Events with SciFi activity
  std::cout << "Done initializing minSciFiHits" << std::endl;
  cutFlow.push_back( new sndAnalysis::sciFiStationCut(0.05, std::vector<int>(1, 1), ch)); // B. Vertex not in 1st wall
  std::cout << "Done initializing SciFi station cut 1" << std::endl;
  cutFlow.push_back( new sndAnalysis::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // C. Vertex not in 5th wall
  std::cout << "Done initializing SciFi station cut 2" << std::endl;
  cutFlow.push_back( new sndAnalysis::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // D. No hits in first SciFi plane
  std::cout << "Done initializing SciFi station cut 3" << std::endl;
  cutFlow.push_back( new sndAnalysis::vetoCut(ch)); // E. No veto hits
  std::cout << "Done initializing veto cut" << std::endl;
  cutFlow.push_back( new sndAnalysis::minSciFiConsecutivePlanes(ch)); // F. At least two consecutive SciFi planes hit
  std::cout << "Done initializing scifi consecutive" << std::endl;
  //  cutFlow.push_back( new sndAnalysis::minSciFiHits(35, ch)); // G. At least 35 hits in SciFi
  cutFlow.push_back( new sndAnalysis::DSActivityCut(ch)); // I. If there is a downstream hit, require hits in all upstream stations.
  std::cout << "Done initializing DS activity" << std::endl;
  //  cutFlow.push_back( new sndAnalysis::USQDCCut(600, ch)); // J. Total US QDC > 600
  if (~isMC) cutFlow.push_back( new sndAnalysis::eventDeltatCut(-1, 100, ch)); // K. Previous event more than 100 clock cycles away. To avoid deadtime issues.
  std::cout << "Done initializing delta t 1" << std::endl;
  if (~isMC) cutFlow.push_back( new sndAnalysis::eventDeltatCut(+1, 10, ch)); // L. Next event more than 10 clock cycles away. To avoid event builder issues.
  std::cout << "Done initializing delta t 2" << std::endl;
  cutFlow.push_back( new sndAnalysis::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // M. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
  std::cout << "Done initializing scifi fid" << std::endl;
  cutFlow.push_back( new sndAnalysis::avgDSFiducialCut(70, 105, 10, 50, ch)); // M. Average SciFi hit channel number must be within [70, 105] (ver) and [10, 50] (hor)
  std::cout << "Done initializing DS fid" << std::endl;

  std::cout << "Done initializing cuts" << std::endl;

  unsigned int n_cuts = cutFlow.size();

  // Book histograms
  // Cut-by-cut
  // All cut variables
  std::vector<std::vector<TH1D*> > cut_by_cut_var_histos = std::vector<std::vector<TH1D*> >();
  for (int i_cut = -1; i_cut < n_cuts; i_cut++){
    std::vector<TH1D*> this_cut_by_cut_var_histos = std::vector<TH1D*>();
    for (sndAnalysis::baseCut * cut : cutFlow) {
      for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
	this_cut_by_cut_var_histos.push_back(new TH1D((std::to_string(i_cut)+"_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
						      cut->getShortName().c_str(), 
						      cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
      }
    }
    cut_by_cut_var_histos.push_back(this_cut_by_cut_var_histos);
  }
  
  std::vector<std::vector<TH1D*> > cut_by_cut_truth_histos = std::vector<std::vector<TH1D*> >();
  if (isMC) {
    for (int i_cut = -1; i_cut < n_cuts; i_cut++){
      std::vector<TH1D*> this_cut_by_cut_truth_histos = std::vector<TH1D*>();
      this_cut_by_cut_truth_histos.push_back(new TH1D((std::to_string(i_cut)+"_Enu").c_str(), "Enu", 300, 0, 3000));
      this_cut_by_cut_truth_histos.push_back(new TH1D((std::to_string(i_cut)+"_EEM").c_str(), "ELep", 300, 0, 3000));
      this_cut_by_cut_truth_histos.push_back(new TH1D((std::to_string(i_cut)+"_EHad").c_str(), "EHad", 300, 0, 3000));
      this_cut_by_cut_truth_histos.push_back(new TH1D((std::to_string(i_cut)+"_vtxX").c_str(), "vtxX", 200, -100, 0));
      this_cut_by_cut_truth_histos.push_back(new TH1D((std::to_string(i_cut)+"_vtxY").c_str(), "vtxY", 200, 0, 100));
      this_cut_by_cut_truth_histos.push_back(new TH1D((std::to_string(i_cut)+"_vtxZ").c_str(), "vtxZ", 200, 280, 380));

      cut_by_cut_truth_histos.push_back(this_cut_by_cut_truth_histos);
    }
  }

  // Effect of a single cut
  std::vector<TH1D*> only_cut_var_histos = std::vector<TH1D*>();
  for (sndAnalysis::baseCut * cut : cutFlow) {
    for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
      only_cut_var_histos.push_back(new TH1D(("only_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
					     cut->getShortName().c_str(), 
					     cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
    }
  }

  // N-1
  std::vector<TH1D*> n_minus_1_var_histos = std::vector<TH1D*>();
  for (sndAnalysis::baseCut * cut : cutFlow) {
    for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
      n_minus_1_var_histos.push_back(new TH1D(("n_minus_1_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
					      cut->getShortName().c_str(), 
					      cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
    }
  }

  std::cout << "Done initializing histograms" << std::endl;

  // Cut flow
  TH1D * cutFlowHistogram = new TH1D("cutFlow", "Cut flow;;Number of events passing cut", n_cuts+1, 0, n_cuts+1);
  for (int i = 2; i <= cutFlowHistogram->GetNbinsX(); i++){
    cutFlowHistogram->GetXaxis()->SetBinLabel(i, cutFlow.at(i-2)->getName().c_str());
  }

  std::cout << "Done initializing cut flow histogram" << std::endl;

  // Get number of entries
  unsigned long int n_entries = ch->GetEntries();

  // Holder for cut results
  std::vector<bool> passes_cut = std::vector<bool>(false, n_cuts);
  int n_cuts_passed = 0;

  std::cout << "Starting event loop" << std::endl;
  for (unsigned long int i_entry = 0; i_entry < n_entries; i_entry++){
    ch->GetEntry(i_entry);
    if (i_entry%10000 == 0) std::cout << "Reading entry " << i_entry << std::endl;
    
    cutFlowHistogram->Fill(0);

    // Apply cuts
    bool accept_event = true;
    int i_cut = 0;
    for (sndAnalysis::baseCut * cut : cutFlow){
      if (cut->passCut()){
	cutFlowHistogram->Fill(1 + i_cut);
	passes_cut[i_cut] = true;
	n_cuts_passed++;
      } else {
	accept_event = false;
	passes_cut[i_cut] = false;
      }
      i_cut++;
    }
    if (accept_event) outTree->Fill();

    // Fill histograms
    std::vector<TH1D*>::iterator hist_it;
    // Sequential
    for (int seq_cut = -1; seq_cut < passes_cut.size(); seq_cut++){
      if (seq_cut >= 0){
	if (passes_cut[seq_cut]) break;
      }
      
      hist_it = cut_by_cut_var_histos[seq_cut+1].begin();
      for (sndAnalysis::baseCut * cut : cutFlow) {
	for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	  (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	  hist_it++;
	}
      }
      cut_by_cut_truth_histos[seq_cut][0]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetEnergy()); // Enu
      cut_by_cut_truth_histos[seq_cut][1]->Fill(((ShipMCTrack*) MCTracks->At(1))->GetEnergy()); // ELep
      cut_by_cut_truth_histos[seq_cut][2]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetEnergy()-((ShipMCTrack*) MCTracks->At(1))->GetEnergy()); // EHad
      cut_by_cut_truth_histos[seq_cut][3]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartX()); // X
      cut_by_cut_truth_histos[seq_cut][4]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartY()); // Y
      cut_by_cut_truth_histos[seq_cut][5]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartZ()); // Z
    }

    // Only cut
    int current_cut = 0;
    hist_it = only_cut_var_histos.begin();
    for (sndAnalysis::baseCut * cut : cutFlow) {
      for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	if (passes_cut[current_cut]) (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	hist_it++;
      }
      current_cut++;
    }
  
    // N-1
    current_cut = 0;
    if (n_cuts_passed == (cutFlow.size() - 1)) {
      hist_it = n_minus_1_var_histos.begin();
      for (sndAnalysis::baseCut * cut : cutFlow) {
	for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	  if (~passes_cut[current_cut]) (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	  hist_it++;
	}
	current_cut++;
      }
    }
  }

  outFile->Write();

  return 0;
}
