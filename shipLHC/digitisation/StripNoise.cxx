#include "StripNoise.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"
#include "TRandom.h"


#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;

StripNoise::StripNoise() {}

AdvSignal StripNoise::AddGaussianNoise(AdvSignal Signal)
{
    std::vector<Int_t> Strips = Signal.getStrips();
    Int_t StripsSize = Strips.size();
    std::vector<Double_t> Amplitude = Signal.getIntegratedSignal();

    TRandom* rndm = gRandom;
    for (int i = 0; i < StripsSize; i++)
    {
        Double_t x = rndm->Gaus(0, stripsensor::frontend::NoiseRMS);
        Amplitude[i] += x ; 
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

    double probabilityLeft = 0.5 * std::erfc(stripsensor::frontend::NoiseSigmaThreshold / std::sqrt(2.0));

    Double_t MeanNumberofNoisyChannels = probabilityLeft * stripsensor::frontend::NumberofStrips; 
    cout << MeanNumberofNoisyChannels << endl ;

    TRandom* rndm = gRandom;
    Double_t NumberofNoisyChannels = rndm->Poisson(MeanNumberofNoisyChannels);

    for (int i = 0; i < NumberofNoisyChannels; i++)
    {
        //need to check if it matches the advsignal strip channels
        NoisyStrips.push_back(rndm->Integer(stripsensor::frontend::NumberofStrips));
        NoisyAmplitudes.push_back(generate_gaussian_tail(stripsensor::frontend::NoiseSigmaThreshold*stripsensor::frontend::NoiseRMS, stripsensor::frontend::NoiseRMS)); 
    }

    for (int j = 0; j < StripsSize; j++)
    {
        Amplitude[j] += rndm->Gaus(0, stripsensor::frontend::NoiseRMS);
    }

    NoiseAddedStrips.insert(NoisyStrips.end(), NoisyStrips.begin(), NoisyStrips.end());
    NoiseAddedAmplitudes.insert(NoisyAmplitudes.end(), NoisyAmplitudes.begin(), NoisyAmplitudes.end());

    AdvSignal NoiseSignal(NoiseAddedStrips, NoiseAddedAmplitudes);

    return NoiseSignal;

}


double StripNoise::generate_gaussian_tail(const double a,const double sigma) 
{
    // taken from CMS 
    
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