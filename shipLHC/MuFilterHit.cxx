#include "MuFilterHit.h"

#include "FairRunSim.h"
#include "MuFilter.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"

#include <TRandom.h>
#include <iomanip>

Double_t speedOfLight = TMath::C() * 100. / 1000000000.0;   // from m/sec to cm/ns
// -----   Default constructor   -------------------------------------------
MuFilterHit::MuFilterHit()
    : SndlhcHit()
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
}
// -----   Standard constructor   ------------------------------------------
MuFilterHit::MuFilterHit(Int_t detID)
    : SndlhcHit(detID)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
}
MuFilterHit::MuFilterHit(Int_t detID, Int_t nP, Int_t nS)
    : SndlhcHit(detID, nP, nS)
{
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
        signals[i] = -999;
        times[i] = -999;
        fDaqID[i] = -1;
    }
}

// -----   constructor from MuFilterPoint   ------------------------------------------
MuFilterHit::MuFilterHit(Int_t detID, std::vector<MuFilterPoint*> V)
    : SndlhcHit()
{
    MuFilter* MuFilterDet = dynamic_cast<MuFilter*>(gROOT->GetListOfGlobals()->FindObject("MuFilter"));
    // get parameters from the MuFilter detector for simulating the digitized information
    nSiPMs = MuFilterDet->GetnSiPMs(detID);
    if (floor(detID / 10000) == 3 && detID % 1000 > 59)
        nSides = MuFilterDet->GetnSides(detID) - 1;
    else
        nSides = MuFilterDet->GetnSides(detID);

    Float_t timeResol = MuFilterDet->GetConfParF("MuFilter/timeResol");

    Float_t attLength = 0;
    Float_t siPMcalibration = 0;
    Float_t siPMcalibrationS = 0;
    Float_t propspeed = 0;
    if (floor(detID / 10000) == 3) {
        if (nSides == 2) {
            attLength = MuFilterDet->GetConfParF("MuFilter/DsAttenuationLength");
        } else {
            attLength = MuFilterDet->GetConfParF("MuFilter/DsTAttenuationLength");
        }
        siPMcalibration = MuFilterDet->GetConfParF("MuFilter/DsSiPMcalibration");
        propspeed = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
    } else {
        attLength = MuFilterDet->GetConfParF("MuFilter/VandUpAttenuationLength");
        siPMcalibration = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibration");
        siPMcalibrationS = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibrationS");
        propspeed = MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");
    }

    for (unsigned int j = 0; j < 16; ++j) {
        signals[j] = -1;
        times[j] = -1;
    }
    LOG(DEBUG) << "detid " << detID << " size " << nSiPMs << "  side " << nSides;

    fDetectorID = detID;
    Float_t signalLeft = 0;
    Float_t signalRight = 0;
    Float_t earliestToAL = 1E20;
    Float_t earliestToAR = 1E20;
    for (auto p = std::begin(V); p != std::end(V); ++p) {

        Double_t signal = (*p)->GetEnergyLoss();

        // Find distances from MCPoint centre to ends of bar
        TVector3 vLeft, vRight;
        TVector3 impact((*p)->GetX(), (*p)->GetY(), (*p)->GetZ());
        MuFilterDet->GetPosition(fDetectorID, vLeft, vRight);
        Double_t distance_Left = (vLeft - impact).Mag();
        Double_t distance_Right = (vRight - impact).Mag();
        signalLeft += signal * TMath::Exp(-distance_Left / attLength);
        signalRight += signal * TMath::Exp(-distance_Right / attLength);

        // for the timing, find earliest particle and smear with time resolution
        Double_t ptime = (*p)->GetTime();
        Double_t t_Left = ptime + distance_Left / propspeed;
        Double_t t_Right = ptime + distance_Right / propspeed;
        if (t_Left < earliestToAL) {
            earliestToAL = t_Left;
        }
        if (t_Right < earliestToAR) {
            earliestToAR = t_Right;
        }
    }
    // shortSiPM = {3,6,11,14,19,22,27,30,35,38,43,46,51,54,59,62,67,70,75,78}; - counting from 1!
    // In the SndlhcHit class the 'signals' array starts from 0.
    for (unsigned int j = 0; j < nSiPMs; ++j) {
        if (j == 2 or j == 5) {
            signals[j] = signalLeft / float(nSiPMs)
                         * siPMcalibrationS;   // most simplest model, divide signal individually. Small SiPMS special
            times[j] = gRandom->Gaus(earliestToAL, timeResol);
        } else {
            signals[j] =
                signalLeft / float(nSiPMs) * siPMcalibration;   // most simplest model, divide signal individually.
            times[j] = gRandom->Gaus(earliestToAL, timeResol);
        }
        if (nSides > 1) {
            signals[j + nSiPMs] =
                signalRight / float(nSiPMs) * siPMcalibration;   // most simplest model, divide signal individually.
            times[j + nSiPMs] = gRandom->Gaus(earliestToAR, timeResol);
        }
    }
    flag = true;
    for (Int_t i = 0; i < 16; i++) {
        fMasked[i] = kFALSE;
    }
    LOG(DEBUG) << "signal created";
}

// -----   Destructor   ----------------------------------------------------
MuFilterHit::~MuFilterHit() {}
// -------------------------------------------------------------------------

// -----   Public method GetEnergy   -------------------------------------------
Float_t MuFilterHit::GetEnergy()
{
    // to be calculated from digis and calibration constants, missing!
    Float_t E = 0;
    for (unsigned int j = 0; j < nSiPMs; ++j) {
        E += signals[j];
        if (nSides > 1) {
            E += signals[j + nSiPMs];
        }
    }
    return E;
}

bool MuFilterHit::isVertical()
{
    if (floor(fDetectorID / 10000) == 3 && fDetectorID % 1000 > 59) {
        return kTRUE;
    } else {
        return kFALSE;
    }
}

bool MuFilterHit::isShort(Int_t i)
{
    if (i % 8 == 2 || i % 8 == 5) {
        return kTRUE;
    } else {
        return kFALSE;
    }
}

// -----   Public method Get List of signals   -------------------------------------------
std::map<Int_t, Float_t> MuFilterHit::GetAllSignals(Bool_t mask, Bool_t positive)
{
    std::map<Int_t, Float_t> allSignals;
    for (unsigned int s = 0; s < nSides; ++s) {
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            unsigned int channel = j + s * nSiPMs;
            if (signals[channel] < -900) {
                continue;
            }
            if (signals[channel] > 0 || !positive) {
                if (!fMasked[channel] || !mask) {
                    allSignals[channel] = signals[channel];
                }
            }
        }
    }
    return allSignals;
}

// -----   Public method Get List of time measurements   -------------------------------------------
std::map<Int_t, Float_t> MuFilterHit::GetAllTimes(Bool_t mask)
{
    std::map<Int_t, Float_t> allTimes;
    for (unsigned int s = 0; s < nSides; ++s) {
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            unsigned int channel = j + s * nSiPMs;
            if (signals[channel] > 0) {
                if (!fMasked[channel] || !mask) {
                    allTimes[channel] = times[channel];
                }
            }
        }
    }
    return allTimes;
}

// -----   Public method Get time difference mean Left - mean Right   -----------------
Float_t MuFilterHit::GetDeltaT(Bool_t mask)
// based on mean TDC measured on Left and Right
{
    Float_t mean[] = {0, 0};
    Int_t count[] = {0, 0};
    Float_t dT = -999.;
    for (unsigned int s = 0; s < nSides; ++s) {
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            unsigned int channel = j + s * nSiPMs;
            if (signals[channel] > 0) {
                if (!fMasked[channel] || !mask) {
                    mean[s] += times[channel];
                    count[s] += 1;
                }
            }
        }
    }
    if (count[0] > 0 && count[1] > 0) {
        dT = mean[0] / count[0] - mean[1] / count[1];
    }
    return dT;
}
Float_t MuFilterHit::GetFastDeltaT(Bool_t mask)
// based on fastest (earliest) TDC measured on Left and Right
{
    Float_t first[] = {1E20, 1E20};
    Float_t dT = -999.;
    for (unsigned int s = 0; s < nSides; ++s) {
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            unsigned int channel = j + s * nSiPMs;
            if (signals[channel] > 0) {
                if (!fMasked[channel] || !mask) {
                    if (times[channel] < first[s]) {
                        first[s] = times[channel];
                    }
                }
            }
        }
    }
    if (first[0] < 1E10 && first[1] < 1E10) {
        dT = first[0] - first[1];
    }
    return dT;
}

// -----   Public method Get mean time  -----------------
Float_t MuFilterHit::GetImpactT(Bool_t mask)
{
    Float_t mean[] = {0, 0};
    Int_t count[] = {0, 0};
    Float_t dT = -999.;
    Float_t dL;
    MuFilter* MuFilterDet = dynamic_cast<MuFilter*>(gROOT->GetListOfGlobals()->FindObject("MuFilter"));
    if (floor(fDetectorID / 10000 == 3)) {
        dL = MuFilterDet->GetConfParF("MuFilterDet/DownstreamBarX") / MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
    } else if (floor(fDetectorID / 10000 == 2)) {
        dL =
            MuFilterDet->GetConfParF("MuFilterDet/UpstreamBarX") / MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");
    } else {
        dL = MuFilterDet->GetConfParF("MuFilterDet/VetoBarX") / MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");
    }

    for (unsigned int s = 0; s < nSides; ++s) {
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            unsigned int channel = j + s * nSiPMs;
            if (signals[channel] > 0) {
                if (!fMasked[channel] || !mask) {
                    mean[s] += times[channel];
                    count[s] += 1;
                }
            }
        }
    }
    if (count[0] > 0 && count[1] > 0) {
        dT = (mean[0] / count[0] + mean[1] / count[1]) / 2. * 6.25 - dL / 2.;   // TDC to ns = 6.25
    }
    return dT;
}

std::map<TString, Float_t> MuFilterHit::SumOfSignals(Bool_t mask)
{
    /*    use cases, for Veto and DS small/large ignored
            sum of signals left large SiPM:    LL
            sum of signals right large SiPM: RL
            sum of signals left small SiPM:    LS
            sum of signals right small SiPM: RS
            sum of signals left and right:
    */
    Float_t theSumL = 0;
    Float_t theSumR = 0;
    Float_t theSumLS = 0;
    Float_t theSumRS = 0;
    for (unsigned int s = 0; s < nSides; ++s) {
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            unsigned int channel = j + s * nSiPMs;
            if (signals[channel] > 0) {
                if (!fMasked[channel] || !mask) {
                    if (s == 0 and !isShort(j)) {
                        theSumL += signals[channel];
                    }
                    if (s == 0 and isShort(j)) {
                        theSumLS += signals[channel];
                    }
                    if (s == 1 and !isShort(j)) {
                        theSumR += signals[channel];
                    }
                    if (s == 1 and isShort(j)) {
                        theSumRS += signals[channel];
                    }
                }
            }
        }
    }
    std::map<TString, Float_t> sumSignals;
    sumSignals["SumL"] = theSumL;
    sumSignals["SumR"] = theSumR;
    sumSignals["SumLS"] = theSumLS;
    sumSignals["SumRS"] = theSumRS;
    sumSignals["Sum"] = theSumL + theSumR;
    sumSignals["SumS"] = theSumLS + theSumRS;
    return sumSignals;
}

// -----   Public method Print   -------------------------------------------
void MuFilterHit::Print() const
{
    std::cout << "-I- MuFilterHit: MuFilter hit " << " in detector " << fDetectorID;

    if (floor(fDetectorID / 10000) == 3 && fDetectorID % 1000 > 59) {
        std::cout << " with vertical bars" << std::endl;
        std::cout << "top digis:";
        for (unsigned int j = 0; j < nSiPMs; ++j) {
            std::cout << signals[j] << " ";
        }
    } else {
        std::cout << " with horizontal bars" << std::endl;
        for (unsigned int s = 0; s < nSides; ++s) {
            if (s == 0) {
                std::cout << "left digis:";
            } else {
                std::cout << "right digis:";
            }
            for (unsigned int j = 0; j < nSiPMs; ++j) {
                std::cout << signals[j] << " ";
            }
        }
    }
    std::cout << std::endl;
}
// -------------------------------------------------------------------------

ClassImp(MuFilterHit)
