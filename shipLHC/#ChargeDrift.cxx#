#include "ChargeDrift.h"
#include "SurfaceSignal.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "SiDigiParameters.h"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std; 

ChargeDrift::ChargeDrift() {}

std::vector<SurfaceSignal> ChargeDrift::Drift(std::vector<EnergyFluctUnit> EnergyLossVector)
{   
    std::vector<SurfaceSignal> DiffusionSignal; 
    for (int k = 0; k < EnergyLossVector.size(); k++)
    {
        vector<Double_t> diffusionarea; 
        vector<TVector3> diffusionpos; 
        vector<Double_t> amplitude; 

        DriftPos = EnergyLossVector[k].getDriftPos();
        EnergyFluct = EnergyLossVector[k].getEfluct();

        for (int i = 0; i < DriftPos.size(); i ++)
        {
            DriftDistance = 0.025 - DriftPos[i].Z();
            DriftDistanceFraction = DriftDistance / stripsensor::drift::module_thickness; 
            DriftDistanceFraction = DriftDistanceFraction > 0. ? DriftDistanceFraction : 0. ; 
            DriftDistanceFraction = DriftDistanceFraction < 1. ? DriftDistanceFraction : 1. ;

            //DriftTime = GetDriftTime(DriftDistance);
            
            Double_t tn = (stripsensor::drift::module_thickness * stripsensor::drift::module_thickness) / (2 * stripsensor::drift::depletion_voltage * stripsensor::drift::charge_mobility);
            DriftTime = -tn * log(1 - 2 * stripsensor::drift::depletion_voltage * DriftDistanceFraction / (stripsensor::drift::depletion_voltage + stripsensor::drift::applied_voltage)) + stripsensor::drift::chargedistributionRMS;

            
            DriftPositiononSurface.SetXYZ(DriftPos[i].X(), DriftPos[i].Y(), DriftPos[i].Z()+DriftDistance); 

            DiffusionConstant = (1.38E-23 / 1.6E-19) * stripsensor::drift::charge_mobility * stripsensor::drift::temperature; 
        
            DiffusionArea = sqrt(2* DiffusionConstant * DriftTime); 
            
            Amplitude = EnergyFluct[i]> 0. ? floor(EnergyFluct[i]*1e9 / stripsensor::drift::perGeV) : 0. ;

            diffusionarea.push_back(DiffusionArea); 
            diffusionpos.push_back(DriftPositiononSurface);
            amplitude.push_back(Amplitude);
            
            // cout << DiffusionArea << endl; 
            // cout << DriftPositiononSurface.X() << endl; 
            
        }

        SurfaceSignal Diffused(diffusionarea, diffusionpos, amplitude);
        
        DiffusionSignal.push_back(Diffused);
    }
    return DiffusionSignal;
}

Double_t ChargeDrift::GetDriftTime(Double_t distance)
{
    Double_t z_width = 400e-6; // in um   
    Double_t p = 120e-6; 
    Double_t w = 30e-6; 
    Double_t h = 5e-6; 

    Double_t Na = 1e15; // in cm 
    Double_t Nd = 1e12; // in cm 

    Double_t e = 1.6e-19; // in C

    Double_t k = 1.38e-23; 
    Double_t T = 300;
    Double_t epsilon = 8.854e-14*11.9; 

    Double_t ni = 1.45e10; // in cm 

    Double_t V_a = 470;
    Double_t V_0 = (((k*T)/e)*log((Na*Nd)/(pow(ni, 2)))) + V_a; 

    Double_t W = sqrt((2*epsilon*(Nd + Na)*V_0)/(e*Nd*Na));

    Double_t Wp = (Nd * W) / (Nd + Na);
    Double_t Wn = (Na * W) / (Nd + Na);

    Int_t num = 80; 

    Double_t x[num] ;
    Double_t y[num] ;
    Double_t E_strip[num] ;

    Double_t x_start = 0 ; // in cm
    Double_t x_end = 500e-4; 


    for (int j = 0; j < num ; j ++ )
    {
        x[j] = x_start + abs(j*((x_start - x_end)/num)); 
        
    }

    for (int i = 0; i < num; i++)
    {
        if ((-Wp <= x[i]) && (x[i] < 0))
        {
            E_strip[i] = -((e*Na*x[i])/epsilon) - ((e*Na*Wp)/epsilon);
            
        }
        else if ((0 <= x[i]) && (x[i] < Wn))
        {
            E_strip[i] = ((e*Nd*x[i])/epsilon) - ((e*Nd*Wn)/epsilon);
        }
        else 
        {
            E_strip[i] = 0;
        }
    }

    Double_t dt[num]; 
    Double_t tfrac[num]; 

    for(int l = 0; l < num ; l++)
    {
        //tfrac[l] = -x[l]/250e-6; 
        dt[l] = -x[l]/(480*E_strip[l]); 
    }

    auto gr = new TGraph (num, x, dt);
    cout << distance << endl; 
    gr->Draw();
    Double_t drifttime = gr->Eval(distance);
    return drifttime; 
}