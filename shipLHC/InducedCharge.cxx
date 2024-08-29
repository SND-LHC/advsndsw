#include "InducedCharge.h"
#include "SurfaceSignal.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TRandom.h"
#include "SiSensor.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std; 

InducedCharge::InducedCharge() {}

void InducedCharge::IntegrateCharge(vector<Int_t> detID, SurfaceSignal DiffusionSignal)
{   
    std::vector<Double_t> amplitude = DiffusionSignal.getAmplitude(); 
    std::vector<TVector3> surfacepos = DiffusionSignal.getSurfacePos();
    std::vector<Double_t> diffusionarea = DiffusionSignal.getDiffusionArea();

    std::vector<Int_t> AffectedStrips; 

    for(int i = 0; i < surfacepos.size(); i++)
    {
        AffectedStrips = GetStrips(surfacepos[i], diffusionarea[i]);
    }
}

std::vector<Int_t> InducedCharge::GetStrips(TVector3 point, Double_t area)
{
    std::vector<Int_t> affectedstrips;

    int fromstrip = floor(((point.X()-(NSigma*area)) / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2));
    fromstrip = max(0, fromstrip);
    fromstrip = min(advsnd::strips - 1, fromstrip);
    affectedstrips.push_back(fromstrip);

    int tostrip = floor(((point.X()+(NSigma*area)) / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2));
    tostrip = max(0, tostrip);
    tostrip = min(advsnd::strips - 1, tostrip);

    affectedstrips.push_back(tostrip);
    return affectedstrips; 
}