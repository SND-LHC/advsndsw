#include "StripNoise.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"
#include "TRandom.h"
#include "TGeoManager.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <bits/stdc++.h>
using namespace std;

/* To Do : 
Add Baseline Shift?
Add Pedestals for each channel after calibration */

StripNoise::StripNoise() {}

std::map<int, int> StripNoise::noise_map;
std::vector<int> StripNoise::detids; 

void StripNoise::CreateDetIDs()
{
    TString path = TString::Format("/cave_1/"
                                   "Detector_0/"
                                   "volAdvTarget_0/");
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    if (nav->CheckPath(path)) {
        nav->cd(path);
        
        TGeoNode *currentnode = nav->GetCurrentNode();

        int nDaughters = currentnode->GetNdaughters();

        for (int n = 0; n < nDaughters; n++)
        {
            TGeoNode* currentdaughternode = currentnode -> GetDaughter(n); 
            string name = currentdaughternode->GetName(); 
            int start_index = name.find("Target_Layer_"); 


            if (start_index == 0){
                int layer = stoi(name.erase(0,13));
                TGeoNode* currentsensornode = currentdaughternode->GetDaughter(0); 
                string sensorname = currentsensornode->GetName(); 
                
                int start_sensor_row_index = sensorname.find("Row_"); 
                int start_sensor_column_index = sensorname.find("Column_"); 
                
                int row = stoi(sensorname.substr(4 + start_sensor_row_index, 1));
                int column = stoi(sensorname.substr(7 + start_sensor_column_index, 1)); 

                int strip = 0; 

                int detid = ((layer )   << 13)
                            | ((row     & 0x3)   << 11)
                            | ((column  & 0x1)   << 10)
                            |  (strip   & 0x3FF);

                detids.push_back(detid);
            } 
            
        }

    } else {
        std::cout << " ====== PATH DOES NOT EXIST  ====== " << path << std::endl; 
    }
}

void StripNoise::CreateNoiseProfile()
{
    ifstream noiseavgfile(stripsensor::frontend::NoiseProfile); 
    if (!noiseavgfile.is_open()) {
        cout << "Noise Profile file doesn't exist!" << endl;
        cout << "Using standard noise from SiDigiParameters.h" << endl; 
        return; 
    } 
    
    std::string line ; 

    while(getline(noiseavgfile, line))

    {
        if (line.empty()) continue; 

        stringstream ss(line);

        string key_str, value_str;

        if ((getline(ss, key_str, ',')) && (getline(ss, value_str, ','))){
            int key = stoi(key_str); 
            int value = stoi(value_str);
            noise_map[key] = value; 
        }
    }

    noiseavgfile.close();
}

AdvSignal StripNoise::AddGaussianNoise(Int_t detID)
{
    std::vector<Int_t> Strips;
    std::vector<Double_t> Amplitude;

    TRandom* rndm = gRandom;
    for (int i = 0; i < stripsensor::frontend::NumberofStrips; i++)
    {
        if (stripsensor::frontend::NoiseOption)
        {
            Double_t x = rndm->Gaus(0, stripsensor::frontend::NoiseRMS);
            Strips.push_back(i); 
            Amplitude.push_back(x);
        }
        else
        {
            Int_t current_detid = (((detID) & ~0x3FF) | i); 
            Strips.push_back(i); 
            if (noise_map[current_detid])
            {
                //Double_t x = rndm->Gaus(0, noise_map[current_detid]); 
                Double_t x = noise_map[current_detid];
                Amplitude.push_back(x);        
            }
            else 
            {
                Amplitude.push_back(0); 
            }
        }
    }
    AdvSignal NoiseSignal(Strips, Amplitude); 
    return NoiseSignal; 
}

AdvSignal StripNoise::AddGaussianTailNoise(AdvSignal Signal)
{
    std::vector<Int_t> Strips = Signal.getStrips();
    Int_t StripsSize = Strips.size();
    std::vector<Double_t> Amplitude = Signal.getIntegratedSignal();

    std::vector<Int_t> NoisyStrips; 
    std::vector<Double_t> NoisyAmplitudes; 

    std::vector<Int_t> NoiseAddedStrips = Strips;  
    std::vector<Double_t> NoiseAddedAmplitudes = Amplitude; 

    TRandom* rndm = gRandom;
    for (int j = 0; j < StripsSize; j++)
    {
        Amplitude[j] += rndm->Gaus(0, stripsensor::frontend::NoiseRMS);
    }

    double probabilityLeft = 0.5 * std::erfc(stripsensor::frontend::NoiseSigmaThreshold / std::sqrt(2.0));

    Double_t MeanNumberofNoisyChannels = probabilityLeft * stripsensor::frontend::NumberofStrips; 
    
    Double_t NumberofNoisyChannels = int(rndm->Poisson(MeanNumberofNoisyChannels)); 
    for (int i = 0; i < NumberofNoisyChannels; i++)
    {
        Strips.push_back(rndm->Integer(stripsensor::frontend::NumberofStrips));
        Amplitude.push_back(generate_gaussian_tail(stripsensor::frontend::NoiseSigmaThreshold*stripsensor::frontend::NoiseRMS, stripsensor::frontend::NoiseRMS)); 
    }
    AdvSignal NoiseSignal(Strips, Amplitude);

    return NoiseSignal;

}

AdvSignal StripNoise::AddCMNoise(AdvSignal Signal)
{
    Double_t CMnoise = stripsensor::frontend::CMNoise;
    std::vector<Int_t> Strips = Signal.getStrips(); 
    std::vector<Double_t> Charge = Signal.getIntegratedSignal(); 
    Int_t nAPVs = int(stripsensor::frontend::NumberofStrips / 128);
    std::vector<Double_t> CMVector ; 
    TRandom* rndm = gRandom;
    for (int i = 0 ; i < stripsensor::frontend::NumberofStrips; i += 128)
    {
        CMVector.push_back(rndm->Gaus(0, CMnoise)); 
    }

    for (int j = 0; j < Strips.size(); j ++)
    {
        Charge[j] += CMVector[int(Strips[j]/128)];
    }

    AdvSignal NoiseSignal(Strips, Charge);
    return NoiseSignal; 
}

double StripNoise::generate_gaussian_tail(const double a,const double sigma) 
{
    // taken from CMSSW 
    
  /* Returns a gaussian random variable larger than a
   * This implementation does one-sided upper-tailed deviates.
   */
    TRandom* rndm = gRandom;
    double s = a / sigma;

    if (s < 1) {
        /*
        For small s, use a direct rejection method. The limit s < 1
        can be adjusted to optimise the overall efficiency
        */
        double x;

        do {
        x = rndm->Gaus(0., 1.0);
        } while (x < s);
        return x * sigma;

    } else {
        /* Use the "supertail" deviates from the last two steps
        * of Marsaglia's rectangle-wedge-tail method, as described
        * in Knuth, v2, 3rd ed, pp 123-128.  (See also exercise 11, p139,
        * and the solution, p586.)
        */

        double u, v, x;

        do {
        u = rndm->Rndm();
        do {
            v = rndm->Rndm();
        } while (v == 0.0);
        x = sqrt(s * s - 2 * log(v));
        } while (x * u > s);
        return x * sigma;
    }
}
