#include "InducedCharge.h"
#include "SurfaceSignal.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TRandom.h"
#include "SiSensor.h"
#include "AdvSignal.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>

using namespace std; 

InducedCharge::InducedCharge() {}

AdvSignal InducedCharge::IntegrateCharge(SurfaceSignal DiffusionSignal)
{   
    std::vector<Double_t> amplitude = DiffusionSignal.getAmplitude(); 
    std::vector<TVector3> surfacepos = DiffusionSignal.getSurfacePos();
    std::vector<Double_t> diffusionarea = DiffusionSignal.getDiffusionArea();

    std::vector<Int_t> AffectedStrips;
    std::vector<Int_t> temp_AffectedStrips; 
    std::vector<Double_t> ChargeDeposited; 
    std::vector<Double_t> TotalChargeDeposited;

    Int_t N ;
    Double_t x_start; 
    Double_t x_end; 
    Double_t z_start; 
    Double_t z_end;  
    Double_t integratedcharge; 

    for (int i = 0; i < surfacepos.size(); i++)
    {
     AffectedStrips = GetStrips(surfacepos[i], diffusionarea[i]);   

         for(int j = 0; j < AffectedStrips.size(); j++)
        {
            x_start = (AffectedStrips[j] - (advsnd::strips / 2))*(advsnd::sensor_width / advsnd::strips) - (strip_pitch / 2); // check calculation
            x_end = (AffectedStrips[j] - (advsnd::strips / 2))*(advsnd::sensor_width / advsnd::strips) + (strip_pitch / 2); // check calculation 
            z_start = (x_start - surfacepos[i].X()) / diffusionarea[j];
            z_end = abs(x_end - surfacepos[i].X()) / diffusionarea[j];
            integratedcharge = (erf((z_end) / TMath::Sqrt2()) / 2) - (erf((z_start) / TMath::Sqrt2()) / 2);
            //cout << surfacepos[j].X() << "\t" << x_start << "\t" << x_end << "\t" << z_start << "\t" << z_end << "\t" << (x_start - surfacepos[j].X()) << "\t" << diffusionarea[j] << "\t" <<integratedcharge << endl; 
            temp_AffectedStrips.push_back(AffectedStrips[j]);
            ChargeDeposited.push_back(integratedcharge*amplitude[j]);        
        }
    }

    std::vector<Int_t> UniqueAffectedStrips = temp_AffectedStrips;  
    sort(UniqueAffectedStrips.begin(), UniqueAffectedStrips.end());
    vector<int>::iterator it;
    it = unique(UniqueAffectedStrips.begin(), UniqueAffectedStrips.end());  

    UniqueAffectedStrips.resize(distance(UniqueAffectedStrips.begin(),it));  
    // cout << "3 : " << AffectedStrips.size() << endl; 
    // cout << "4 : " << surfacepos.size() << endl; 
    // for(int j = 0; j < AffectedStrips.size(); j++)
    // {
    //     x_start = (AffectedStrips[j] - (advsnd::strips / 2))*(advsnd::sensor_width / advsnd::strips) - (strip_pitch / 2); // check calculation
    //     x_end = (AffectedStrips[j] - (advsnd::strips / 2))*(advsnd::sensor_width / advsnd::strips) + (strip_pitch / 2); // check calculation 
    //     z_start = (x_start - surfacepos[j].X()) / diffusionarea[j];
    //     z_end = abs(x_end - surfacepos[j].X()) / diffusionarea[j];
    //     integratedcharge = (erf((z_end) / TMath::Sqrt2()) / 2) - (erf((z_start) / TMath::Sqrt2()) / 2);
    //     //cout << surfacepos[j].X() << "\t" << x_start << "\t" << x_end << "\t" << z_start << "\t" << z_end << "\t" << (x_start - surfacepos[j].X()) << "\t" << diffusionarea[j] << "\t" <<integratedcharge << endl; 
    //     ChargeDeposited.push_back(integratedcharge*amplitude[j]);        
    // }

    Double_t z = accumulate(amplitude.begin(), amplitude.end(), 0);

    for (int k = 0; k < UniqueAffectedStrips.size(); k++)
    {
        Double_t temp_totalcharge = 0.0;
        auto it = find(temp_AffectedStrips.begin(), temp_AffectedStrips.end(), UniqueAffectedStrips[k]); 
        while (it != temp_AffectedStrips.end()) { 
            temp_totalcharge = temp_totalcharge + ChargeDeposited[it - temp_AffectedStrips.begin()];
            it = find(it + 1, temp_AffectedStrips.end(), UniqueAffectedStrips[k]); 
        } 
        TotalChargeDeposited.push_back(temp_totalcharge);
    }

    Double_t r = accumulate(TotalChargeDeposited.begin(), TotalChargeDeposited.end(), 0);

    Double_t rescale_ratio = r/z; 
    for (int k = 0; k < UniqueAffectedStrips.size(); k++)
    {
        TotalChargeDeposited[k] = TotalChargeDeposited[k] * rescale_ratio ; 
    }

    AdvSignal ResponseSignal(UniqueAffectedStrips, TotalChargeDeposited); 
    return ResponseSignal; 
}

std::vector<Int_t> InducedCharge::GetStrips(TVector3 point, Double_t area)
{
    std::vector<Int_t> affectedstrips;

    // for(int i = 0; i < point.size(); i++)
    // {
        int fromstrip = floor(((point.X()-(NSigma*area)) / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2)); // check calculation 
        fromstrip = max(0, fromstrip);
        fromstrip = min(advsnd::strips - 1, fromstrip); 

        //cout << area[i] << "\t" << point[i].X()-(NSigma*area[i]) << endl; 

        int tostrip = floor(((point.X()+(NSigma*area)) / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2));
        tostrip = max(0, tostrip);
        tostrip = min(advsnd::strips - 1, tostrip);

        Int_t N; 
        N = tostrip - fromstrip; 

        //cout << N << "\t" << fromstrip << "\t" << tostrip << "\t" << point[i].X()-(NSigma*area[i]) << "\t" << point[i].X()+(NSigma*area[i]) << endl; 

        for (int i = 0 ; i <= N ; i++)
        {
            affectedstrips.push_back(fromstrip + i);
        }
    // }
    return affectedstrips; 
}