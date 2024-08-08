#include "AdvTargetHit.h"
#include "EnergyFluctUnit.h"

#include "AdvTargetPoint.h"
#include "ChargeDivision.h"
#include "FairLogger.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"
#include "TRandom.h"

#include <TDatabasePDG.h>
#include <iomanip>
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
    ChargeDivision chargedivision{};
    std::string inputfile =
        "advsndsw/shipLHC/data/APVShapeDeco_default.txt";   // change this full path in configuration file
    chargedivision.ReadPulseShape(inputfile);
    EnergyFluctUnit EnergyLossVector = chargedivision.Divide(detID, V);
    EFluct = EnergyLossVector.getEfluct();
    segLen = EnergyLossVector.getsegLen();
    if (EFluct.empty()) {
        EFluctSize = 0;
    } else {
        EFluctSize = size(EFluct);
    }

    // std::string inputfile =
    //     "advsndsw/shipLHC/data/APVShapeDeco_default.txt";   // change this full path in configuration file
    // chargedivision.ReadPulseShape(inputfile);

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
