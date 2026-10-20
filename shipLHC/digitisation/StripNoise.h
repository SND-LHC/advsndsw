#ifndef SHIPLHC_STRIPNOISE_H_
#define SHIPLHC_STRIPNOISE_H_

#include "StripNoise.h"
#include "AdvSignal.h"

#include <iostream>
#include <vector>
#include <map>

class StripNoise
{
  public:
    StripNoise();

    static std::vector<int> detids; 
    static std::map<int, int> noise_map;
    static std::map<int, std::vector<int>> pedestal_map; 

    void CreateDetIDs();
    void CreateNoiseProfile(); 
    void CreatePedestalProfile();
    
    AdvSignal AddGaussianNoise(Int_t detID);
    AdvSignal AddPedestals(Int_t detID);
    AdvSignal AddGaussianTailNoise(AdvSignal Signal);
    AdvSignal AddCMNoise(AdvSignal Signal);
    double generate_gaussian_tail(const double a, const double sigma); 

    void TestingNoise();

};

#endif  