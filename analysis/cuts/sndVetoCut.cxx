#include "sndVetoCut.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"

namespace sndAnalysis {

  vetoCut::vetoCut(TChain * ch) : MuFilterBaseCut(ch) {
    cutName = "No hits in veto";

    shortName = "NoVetoHits";
    nbins = std::vector<int>{16};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{16};
    plot_var = std::vector<double>{-1};

  }

  bool vetoCut::passCut(){
    MuFilterHit * hit;
    TIter hitIterator(muFilterDigiHitCollection);
    
    plot_var[0] = 0;


    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->GetSystem() == 1) plot_var[0] += 1;
    }

    if (plot_var[0] > 0) return false;
    return true;
  }
}
