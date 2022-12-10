#include "sndMuFilterBaseCut.h"

#include <vector>

#include "TClonesArray.h"
#include "MuFilterHit.h"

namespace sndAnalysis {

  TClonesArray * MuFilterBaseCut::muFilterDigiHitCollection = 0;

  MuFilterBaseCut::MuFilterBaseCut(TChain * ch){
    if (muFilterDigiHitCollection == 0){
      muFilterDigiHitCollection = new TClonesArray("MuFilterHit", 470);
      ch->SetBranchAddress("Digi_MuFilterHits", &muFilterDigiHitCollection); 
    }
  }
}
