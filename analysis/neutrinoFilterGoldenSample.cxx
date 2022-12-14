#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"

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
  TChain * ch = new TChain("rawConv");
  ch->Add(argv[1]);
  std::cout << "Got input tree" << std::endl;

  // Output file
  TFile * outFile = new TFile(argv[2], "RECREATE");
  std::cout << "Got output file" << std::endl;  

  // Set up all branches to copy to output TTree.
  TTree * outTree = ch->CloneTree(0);
  std::cout << "Got output tree" << std::endl;  

  // Book histograms

  // Three sets of histograms: only this cut; N-1; sequential. One histogram per cut.

  // Set up cuts
  std::vector< sndAnalysis::baseCut * > cutFlow;
  cutFlow.push_back( new sndAnalysis::minSciFiHits(1, ch)); // A. Events with SciFi activity
  cutFlow.push_back( new sndAnalysis::sciFiStationCut(0.05, std::vector<int>(1, 1), ch)); // B. Vertex not in 1st wall
  cutFlow.push_back( new sndAnalysis::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // C. Vertex not in 5th wall
  cutFlow.push_back( new sndAnalysis::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // D. No hits in first SciFi plane
  cutFlow.push_back( new sndAnalysis::vetoCut(ch)); // E. No veto hits
  cutFlow.push_back( new sndAnalysis::minSciFiConsecutivePlanes(ch)); // F. At least two consecutive SciFi planes hit
  cutFlow.push_back( new sndAnalysis::minSciFiHits(35, ch)); // G. At least 35 hits in SciFi
  cutFlow.push_back( new sndAnalysis::DSActivityCut(ch)); // I. If there is a downstream hit, require hits in all upstream stations.
  cutFlow.push_back( new sndAnalysis::USQDCCut(600, ch)); // J. Total US QDC > 600
  cutFlow.push_back( new sndAnalysis::eventDeltatCut(-1, 100, ch)); // K. Previous event more than 100 clock cycles away. To avoid deadtime issues.
  cutFlow.push_back( new sndAnalysis::eventDeltatCut(+1, 10, ch)); // L. Next event more than 10 clock cycles away. To avoid event builder issues.
  cutFlow.push_back( new sndAnalysis::avgSciFiFiducialCut(200, 1200, 300, 128*12, ch)); // M. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max] (hor)
  cutFlow.push_back( new sndAnalysis::avgDSFiducialCut(70, 105, 10, 50, ch)); // M. Average SciFi hit channel number must be within [70, 110] (ver) and [10, 50] (hor)

  unsigned int n_cuts = cutFlow.size();

  // Cut flow
  TH1D * cutFlowHistogram = new TH1D("cutFlow", "Cut flow;;Number of events passing cut", n_cuts+1, 0, n_cuts+1);
  for (int i = 2; i <= cutFlowHistogram->GetNbinsX(); i++){
    cutFlowHistogram->GetXaxis()->SetBinLabel(i, cutFlow.at(i-2)->getName().c_str());
  }

  // Get number of entries
  unsigned long int n_entries = ch->GetEntries();

  // Holder for cut results
  std::vector<bool> passes_cut = std::vector<bool>(false, n_cuts);

  for (unsigned long int i_entry = 0; i_entry < n_entries; i_entry++){
    ch->GetEntry(i_entry);
    
    cutFlowHistogram->Fill(0);

    //    int n_pass = 0;
    bool accept_event = true;
    int i_cut = 0;
    for (sndAnalysis::baseCut * cut : cutFlow){
      if (cut->passCut()){
	cutFlowHistogram->Fill(1 + i_cut++);
      } else {
	accept_event = false;
	break;
      }
    }
    if (accept_event) outTree->Fill();


  }

  outFile->Write();

  return 0;
}
