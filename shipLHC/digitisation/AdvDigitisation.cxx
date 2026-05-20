#include "AdvDigitisation.h"
#include "AdvTargetPoint.h"
#include "ChargeDivision.h"
#include "ChargeDrift.h"
#include "InducedCharge.h"
#include "FrontendDriver.h"
#include "StripNoise.h"
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
#include <map>
#include <algorithm>
using namespace std;

// Running the digitisation 

/* To be included :
    Saturation of FED dynamic range depending on mode 
    Add noise
    Add FED modes */

AdvDigitisation::AdvDigitisation() {}

std::map<std::string, std::vector<Int_t>> AdvDigitisation::digirunoutput(Int_t detID, const std::vector<AdvTargetPoint *> &V)
{
    // Charge Division
    ChargeDivision chargedivision{};
    std::vector<EnergyFluctUnit> EnergyLossVector = chargedivision.Divide(detID, V);

    //Charge Drift
    ChargeDrift chargedrift{};
    std::vector<SurfaceSignal> DiffusionSignal = chargedrift.Drift(EnergyLossVector);

    //Induced Charge on strips
    InducedCharge inducedcharge{};
    AdvSignal ResponseSignal = inducedcharge.IntegrateCharge(DiffusionSignal);

    //FED Response 
    FrontendDriver frontenddriver{};
    AdvSignal FEDResponseSignal = frontenddriver.FEDResponse(ResponseSignal);

    //Creating map of hit 
    std::map<std::string, std::vector<Int_t>> DigitisedHit; 
    std::vector<Double_t> Charge = FEDResponseSignal.getIntegratedSignal();
    std::vector<Int_t> Strips = FEDResponseSignal.getStrips();
    std::vector<Int_t> ADC(Charge.size()); 
    std::transform(Charge.begin(), Charge.end(), ADC.begin(), [](Double_t x) { 
        if (x>0){return (int)x;}
        else {return 0;}
        });
    DigitisedHit["Strips"] = Strips; 
    DigitisedHit["ADC"] = ADC; 
    
    return DigitisedHit;
}