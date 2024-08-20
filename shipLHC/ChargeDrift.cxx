#include "ChargeDrift.h"
#include "SurfaceSignal.h"
#include "TVector3.h"

#include <iostream>
#include <vector>

using namespace std; 

ChargeDrift::ChargeDrift() {}

SurfaceSignal ChargeDrift::Drift(EnergyFluctUnit EnergyLossVector)
{   
    vector<Double_t> diffusionarea; 
    vector<TVector3> diffusionpos; 

    DriftPos = EnergyLossVector.getDriftPos();
    for (int i = 0; i < DriftPos.size(); i ++)
    {
        DriftDistance = 0 - DriftPos[i].Z();
        DriftDistanceFraction = DriftDistance / module_thickness; 
        DriftDistanceFraction = DriftDistanceFraction > 0. ? DriftDistanceFraction : 0. ; 
        DriftDistanceFraction = DriftDistanceFraction < 1. ? DriftDistanceFraction : 1. ;

        Double_t tn = (module_thickness * module_thickness) / (2 * depletion_voltage * charge_mobility);
        DriftTime = -tn * log(1 - 2 * depletion_voltage * DriftDistanceFraction / (depletion_voltage + applied_voltage)) + chargedistributionRMS;

        DriftPositiononSurface.SetXYZ(DriftPos[i].X(), DriftPos[i].Y(), DriftPos[i].Z()+DriftDistance);

        DiffusionConstant = (1.38E-23 / 1.6E-19) * charge_mobility * temperature; 
       
        DiffusionArea = sqrt(2* DiffusionConstant * DriftTime); 

        diffusionarea.push_back(DiffusionArea); 
        diffusionpos.push_back(DriftPositiononSurface);
        
        // cout << DiffusionArea << endl; 
        // cout << DriftPositiononSurface.X() << endl; 
        
    }

    SurfaceSignal DiffusionSignal(diffusionarea, diffusionpos);
    
    return DiffusionSignal; 
    
}