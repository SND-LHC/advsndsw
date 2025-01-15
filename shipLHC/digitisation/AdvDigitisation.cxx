#include "AdvDigitisation.h"
#include "AdvTargetPoint.h"
#include "ChargeDivision.h"
#include "ChargeDrift.h"
#include "InducedCharge.h"
#include "FrontendDriver.h"
#include "Clustering.h"
#include "EnergyFluctUnit.h"
#include "SurfaceSignal.h"
#include "AdvSignal.h"
#include "SiDigiParameters.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TStopwatch.h"

#include <TSystem.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

// Running the digitisation 

AdvDigitisation::AdvDigitisation() {}

std::vector<AdvSignal> AdvDigitisation::digirunoutput(Int_t detID, const std::vector<AdvTargetPoint *> &V)
{
    // Charge Division
    ChargeDivision chargedivision{};
    std::vector<EnergyFluctUnit> EnergyLossVector = chargedivision.Divide(detID, V);

    //Charge Drift
    ChargeDrift chargedrift{};
    std::vector<SurfaceSignal> DiffusionSignal = chargedrift.Drift(EnergyLossVector);

    //Induced Charge on strips
    InducedCharge inducedcharge{};
    std::vector<AdvSignal> ResponseSignal = inducedcharge.IntegrateCharge(DiffusionSignal);

    return ResponseSignal;
}