#include "sndDSActivityCut.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"

#include <vector>
#include <numeric>

namespace sndAnalysis {

  DSActivityCut::DSActivityCut(TChain * ch) : MuFilterBaseCut(ch) {cutName = "If there are DS hits, all US planes must be hit";}

  bool DSActivityCut::passCut(){
    MuFilterHit * hit;
    TIter hitIterator(muFilterDigiHitCollection);

    bool ds = false;
    std::vector<bool> us = std::vector<bool>(5, false); 

    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->GetSystem() == 2) {
	us[hit->GetPlane()] = true;
      } else if (hit->GetSystem() == 3) {
	ds = true;
      }
    }
    
    if (ds and (std::accumulate(us.begin(), us.end(), 0) == 5)) return true;
    return false;
  }
}
