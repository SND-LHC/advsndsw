#ifndef SHIPLHC_INDUCEDCHARGE_H_
#define SHIPLHC_INDUCEDCHARGE_H_

#include "EnergyFluctUnit.h"
#include "SurfaceSignal.h"
#include "AdvSignal.h"

#include <iostream>
#include <vector>

class InducedCharge
{
  public:
    InducedCharge();
    std::vector<AdvSignal> IntegrateCharge(std::vector<SurfaceSignal> DiffusionSignal);
    std::vector<Int_t> GetStrips(TVector3 point, Double_t area);
    std::vector<std::vector<Double_t>> GetPulseShape(std::string PulseFileName, std::vector<Double_t> ChargeDeposited);
    AdvSignal Coupling(std::vector<Double_t> TotalCharge, std::vector<Int_t> AffectedStrips, std::vector<std::vector<Double_t>> PulseResponse);

};

#endif  