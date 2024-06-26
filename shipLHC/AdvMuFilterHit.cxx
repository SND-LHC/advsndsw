#include "AdvMuFilterHit.h"

#include "FairLogger.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"

#include <TRandom.h>
#include <iomanip>

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
    // MuFilter* MuFilterDet = dynamic_cast<MuFilter*>(gROOT->GetListOfGlobals()->FindObject("MuFilter"));
    // // get parameters from the MuFilter detector for simulating the digitized information
    // nSiPMs = MuFilterDet->GetnSiPMs(detID);
    // if (floor(detID / 10000) == 3 && detID % 1000 > 59)
    //     nSides = MuFilterDet->GetnSides(detID) - 1;
    // else
    //     nSides = MuFilterDet->GetnSides(detID);

    // Float_t timeResol = MuFilterDet->GetConfParF("MuFilter/timeResol");

    // Float_t attLength = 0;
    // Float_t siPMcalibration = 0;
    // Float_t siPMcalibrationS = 0;
    // Float_t propspeed = 0;
    // if (floor(detID / 10000) == 3) {
    //     if (nSides == 2) {
    //         attLength = MuFilterDet->GetConfParF("MuFilter/DsAttenuationLength");
    //     } else {
    //         attLength = MuFilterDet->GetConfParF("MuFilter/DsTAttenuationLength");
    //     }
    //     siPMcalibration = MuFilterDet->GetConfParF("MuFilter/DsSiPMcalibration");
    //     propspeed = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
    // } else {
    //     attLength = MuFilterDet->GetConfParF("MuFilter/VandUpAttenuationLength");
    //     siPMcalibration = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibration");
    //     siPMcalibrationS = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibrationS");
    //     propspeed = MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");
    // }

    // for (unsigned int j = 0; j < 16; ++j) {
    //     signals[j] = -1;
    //     times[j] = -1;
    // }
    // LOG(DEBUG) << "detid " << detID << " size " << nSiPMs << "  side " << nSides;

    // fDetectorID = detID;
    // Float_t signalLeft = 0;
    // Float_t signalRight = 0;
    // Float_t earliestToAL = 1E20;
    // Float_t earliestToAR = 1E20;
    // for (auto p = std::begin(V); p != std::end(V); ++p) {

    //     Double_t signal = (*p)->GetEnergyLoss();

    //     // Find distances from MCPoint centre to ends of bar
    //     TVector3 vLeft, vRight;
    //     TVector3 impact((*p)->GetX(), (*p)->GetY(), (*p)->GetZ());
    //     MuFilterDet->GetPosition(fDetectorID, vLeft, vRight);
    //     Double_t distance_Left = (vLeft - impact).Mag();
    //     Double_t distance_Right = (vRight - impact).Mag();
    //     signalLeft += signal * TMath::Exp(-distance_Left / attLength);
    //     signalRight += signal * TMath::Exp(-distance_Right / attLength);

    //     // for the timing, find earliest particle and smear with time resolution
    //     Double_t ptime = (*p)->GetTime();
    //     Double_t t_Left = ptime + distance_Left / propspeed;
    //     Double_t t_Right = ptime + distance_Right / propspeed;
    //     if (t_Left < earliestToAL) {
    //         earliestToAL = t_Left;
    //     }
    //     if (t_Right < earliestToAR) {
    //         earliestToAR = t_Right;
    //     }
    // }
    // // shortSiPM = {3,6,11,14,19,22,27,30,35,38,43,46,51,54,59,62,67,70,75,78}; - counting from 1!
    // // In the SndlhcHit class the 'signals' array starts from 0.
    // for (unsigned int j = 0; j < nSiPMs; ++j) {
    //     if (j == 2 or j == 5) {
    //         signals[j] = signalLeft / float(nSiPMs)
    //                      * siPMcalibrationS;   // most simplest model, divide signal individually. Small SiPMS
    //                      special
    //         times[j] = gRandom->Gaus(earliestToAL, timeResol);
    //     } else {
    //         signals[j] =
    //             signalLeft / float(nSiPMs) * siPMcalibration;   // most simplest model, divide signal individually.
    //         times[j] = gRandom->Gaus(earliestToAL, timeResol);
    //     }
    //     if (nSides > 1) {
    //         signals[j + nSiPMs] =
    //             signalRight / float(nSiPMs) * siPMcalibration;   // most simplest model, divide signal individually.
    //         times[j + nSiPMs] = gRandom->Gaus(earliestToAR, timeResol);
    //     }
    // }
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
