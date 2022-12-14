#include "sndVetoCut.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"

namespace sndAnalysis {

  vetoCut::vetoCut(TChain * ch) : MuFilterBaseCut(ch) {cutName = "No hits in veto";}

  bool vetoCut::passCut(){
    MuFilterHit * hit;
    TIter hitIterator(muFilterDigiHitCollection);
    
    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->GetSystem() == 1) return false;
    }
    return true;
  }
}
