#include "sndMuFilterBaseCut.h"

#include "MuFilterHit.h"
#include "TClonesArray.h"

#include <vector>

namespace sndAnalysis {

TClonesArray* MuFilterBaseCut::muFilterDigiHitCollection = 0;

MuFilterBaseCut::MuFilterBaseCut(TChain* ch)
{
    if (muFilterDigiHitCollection == 0) {
        muFilterDigiHitCollection = new TClonesArray("MuFilterHit", 470);
        ch->SetBranchAddress("Digi_MuFilterHits", &muFilterDigiHitCollection);
    }
}
}   // namespace sndAnalysis
