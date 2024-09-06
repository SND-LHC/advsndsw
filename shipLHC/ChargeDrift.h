#ifndef SHIPLHC_CHARGEDRIFT_H_
#define SHIPLHC_CHARGEDRIFT_H_

#include "EnergyFluctUnit.h"
#include "SurfaceSignal.h"

#include <iostream>
#include <vector>

class ChargeDrift
{
  public:
    ChargeDrift();
    std::vector<SurfaceSignal> Drift(std::vector<EnergyFluctUnit> EnergyLossVector);
    Double_t GetDriftTime(Double_t distance);

  protected:
    Double_t module_thickness = 0.05; //check value
    int depletion_voltage = 170; 
    int applied_voltage = 300; 
    int charge_mobility = 480;
    Double_t chargedistributionRMS = 6.5e-10;
    int temperature = 300; 
    Double_t perGeV = 3.61; 


    std::vector<TVector3> DriftPos; 
    std::vector<Double_t> EnergyFluct; 
    Double_t DriftDistance; 
    Double_t DriftDistanceFraction;
    Double_t DriftTime; 
    TVector3 DriftPositiononSurface; 
    Double_t DiffusionArea;
    Double_t DiffusionConstant;
    Double_t Amplitude; 

};

#endif   // SHIPLHC_CHARGEDIVISION_H_