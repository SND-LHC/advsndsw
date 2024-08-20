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
    SurfaceSignal Drift(EnergyFluctUnit EnergyLossVector);

  protected:
    Double_t module_thickness = 0.025; //check value
    int depletion_voltage = 170; 
    int applied_voltage = 300; 
    int charge_mobility = 310;
    Double_t chargedistributionRMS = 6.5e-10;
    int temperature = 273; 


    std::vector<TVector3> DriftPos; 
    Double_t DriftDistance; 
    Double_t DriftDistanceFraction;
    Double_t DriftTime; 
    TVector3 DriftPositiononSurface; 
    Double_t DiffusionArea;
    Double_t DiffusionConstant;

};

#endif   // SHIPLHC_CHARGEDIVISION_H_