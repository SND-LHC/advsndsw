#include "sndScifiHit.h"

#include "FairRunSim.h"
#include "Scifi.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"

#include <TRandom.h>
#include <iomanip>

Float_t sndScifiHit::MeanAndRMS(Double_t ly, Float_t nphe_max)
{
    Double_t mean = 0.274572 + 0.980999 * ly - 0.00435048 * ly * ly + 1.0901e-05 * pow(ly, 3) - 1.47437e-08 * pow(ly, 4)
                    + 8.28684e-12 * pow(ly, 5);
    Double_t RMS = -0.00055223 + 0.0667554 * ly - 0.000478245 * ly * ly + 1.48809e-06 * pow(ly, 3)
                   - 2.26064e-09 * pow(ly, 4) + 1.35786e-12 * pow(ly, 5);
    Float_t r = gRandom->Gaus(mean, RMS);
    Float_t signal = std::min(nphe_max, r);
    return signal;
}

// -----   Default constructor   -------------------------------------------
sndScifiHit::sndScifiHit()
    : SndlhcHit()
{
    flag = true;
}
// -----   Standard constructor   ------------------------------------------
sndScifiHit::sndScifiHit(Int_t detID)
    : SndlhcHit(detID)
{
    flag = true;
}
// -----   constructor from point class  ------------------------------------------
sndScifiHit::sndScifiHit(int SiPMChan, std::vector<ScifiPoint*> V, std::vector<Float_t> W)
{
    Scifi* ScifiDet = dynamic_cast<Scifi*>(gROOT->GetListOfGlobals()->FindObject("Scifi"));
    Float_t nphe_min = ScifiDet->GetConfParF("Scifi/nphe_min");
    Float_t nphe_max = ScifiDet->GetConfParF("Scifi/nphe_max");
    Float_t timeResol = ScifiDet->GetConfParF("Scifi/timeResol");
    Float_t signalSpeed = ScifiDet->GetConfParF("Scifi/signalSpeed");

    nSides = 1;
    nSiPMs = 1;
    for (unsigned int j = 0; j < 16; ++j) {
        signals[j] = -1;
        times[j] = -1;
    }

    fDetectorID = SiPMChan;
    Float_t signalTotal = 0;
    Float_t earliestToA = 1E20;
    for (int i = 0; i < V.size(); i++) {

        Double_t signal = V[i]->GetEnergyLoss() * W[i];

        // Find distances from MCPoint centre to ends of fibre
        TVector3 a, b;
        TVector3 impact(V[i]->GetX(), V[i]->GetY(), V[i]->GetZ());
        ScifiDet->GetSiPMPosition(V[i]->GetDetectorID(), a, b);

        Double_t distance;

        bool verticalHit = int(fDetectorID / 100000) % 10 == 1;

        // Calculate distance from energy deposit to SiPM.
        // Vertical
        if (verticalHit) {
            distance = (b - impact).Mag();
            // Horizontal
        } else {
            distance = (impact - a).Mag();
        }

        signalTotal += signal * ly_loss_mean(distance, ly_loss_params);

        // for the timing, find earliest light to arrive at SiPM and smear with time resolution
        Float_t arrival_time = V[i]->GetTime() + distance / signalSpeed;
        if (arrival_time < earliestToA) {
            earliestToA = arrival_time;
        }
    }
    Float_t ly = signalTotal / 0.180 * 1000. * 20.;   // 20 p.e. per 180 keV
    signals[0] = MeanAndRMS(ly, nphe_max);
    if (ly > nphe_min) {   // nominal threshold at 3.5 p.e.
        flag = true;
    } else {
        flag = false;
    }
    times[0] = gRandom->Gaus(earliestToA, timeResol);

    LOG(DEBUG) << "signal created";
}

// -----   Destructor   ----------------------------------------------------
sndScifiHit::~sndScifiHit() {}
// -------------------------------------------------------------------------

// -----   Public method GetEnergy   -------------------------------------------
Float_t sndScifiHit::GetEnergy()
{
    // to be calculated from digis and calibration constants, missing!
    return signals[0];
}

Float_t sndScifiHit::ly_loss_mean(Float_t distance, Float_t* params)
{
    //	It returns the light yield depending on the distance to SiPM
    Double_t att =
        params[0] * TMath::Exp(params[1] * distance / 100.) + params[2] * TMath::Exp(params[3] * distance / 100.);
    return att / (params[0] + params[2]);
}

// -----   Public method Print   -------------------------------------------
void sndScifiHit::Print()
{
    std::cout << "-I- sndScifiHit: Scifi hit " << " in station " << GetStation();
    if (isVertical()) {
        std::cout << " vertical plane ";
    } else {
        std::cout << " horizontal plane ";
    }
    std::cout << "SiPM nr " << GetSiPM() << " SiPM channel " << GetSiPMChan() << std::endl;
}
// -------------------------------------------------------------------------

ClassImp(sndScifiHit)
