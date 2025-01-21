#include "AdvTargetHit.h"
#include "AdvTargetPoint.h"
#include "digitisation/AdvSignal.h"
#include "digitisation/AdvDigitisation.h"
#include "FairLogger.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TStopwatch.h"

#include <TDatabasePDG.h>
#include <iomanip>
#include <typeinfo>
#include <iostream>
#include <string>
using std::cout;
using std::endl;

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
    AdvDigitisation advdigi{};
    AdvSignal DigitisedHit = advdigi.digirunoutput(detID, V);
    fStrips = DigitisedHit.getStrips();
    fCharge = DigitisedHit.getIntegratedSignal();
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
