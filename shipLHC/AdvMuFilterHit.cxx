#include "AdvMuFilterHit.h"

#include "FairLogger.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"

// -----   Default constructor   -------------------------------------------
AdvMuFilterHit::AdvMuFilterHit()
    : SndlhcHit()
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
}
// -----   Standard constructor   ------------------------------------------
AdvMuFilterHit::AdvMuFilterHit(Int_t detID)
    : SndlhcHit(detID)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
}

// -----   constructor from AdvMuFilterPoint   ------------------------------------------
AdvMuFilterHit::AdvMuFilterHit(Int_t detID, const std::vector<AdvMuFilterPoint*>& V)
    : SndlhcHit(detID)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
    LOG(DEBUG) << "signal created";
}

// -----   Public method Print   -------------------------------------------
void AdvMuFilterHit::Print() const
{
    LOG(INFO) << " AdvMuFilterHit: AdvMuFilter hit "
              << " in detector " << fDetectorID;
}
