#include "AdvTargetHit.h"

#include "FairLogger.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"

// -----   Default constructor   -------------------------------------------
AdvTargetHit::AdvTargetHit()
    : SndlhcHit()
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
}
// -----   Standard constructor   ------------------------------------------
AdvTargetHit::AdvTargetHit(Int_t detID)
    : SndlhcHit(detID)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
}

// -----   constructor from AdvMuFilterPoint   ------------------------------------------
AdvTargetHit::AdvTargetHit(Int_t detID, const std::vector<AdvTargetPoint*>& V)
    : SndlhcHit(detID)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
    LOG(DEBUG) << "signal created";
}

AdvTargetHit::AdvTargetHit(Int_t detID, const std::vector<const AdvTargetPoint*>& V)
    : SndlhcHit(detID)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
    LOG(DEBUG) << "signal created";
}

// -----   Public method Print   -------------------------------------------
void AdvTargetHit::Print() const
{
    LOG(INFO) << " AdvTargetHit: AdvMuFilter hit "
              << " in detector " << fDetectorID;
}
