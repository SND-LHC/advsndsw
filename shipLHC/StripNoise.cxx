#include "StripNoise.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"
#include "TRandom.h"

#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_result.h>

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

    gsl_sf_result result;
    int status = gsl_sf_erf_Q_e(stripsensor::frontend::NoiseSigmaThreshold, &result);

    float probabilityLeft = result.val;

    Double_t MeanNumberofNoisyChannels = probabilityLeft * StripsSize; 

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
  /* Returns a gaussian random variable larger than a
   * This implementation does one-sided upper-tailed deviates.
   */
    TRandom* rndm = gRandom;
    double s = a / sigma;
    cout << s << endl; 

    if (s < 1) {
        /*
        For small s, use a direct rejection method. The limit s < 1
        can be adjusted to optimise the overall efficiency
        */
        double x;

        do {
        x = rndm->Gaus(0., 1.0);
        cout << x << endl; 
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
        cout << u << endl; 
        do {
            v = rndm->Rndm();
            cout << v << endl; 
        } while (v == 0.0);
        x = sqrt(s * s - 2 * log(v));
        cout << "2 : " << x << endl; 
        } while (x * u > s);
        return x * sigma;
    }
}

void TestingGaussianNoise()
{
    std::vector<Int_t> Strips = {100, 101, 102}; 
    std::vector<Double_t> IntegratedSignal = {200, 600, 50};
    
    AdvSignal TestSignal(Strips, IntegratedSignal); 

    StripNoise stripnoise; 
    // stripnoise.AddGaussianNoise(TestSignal);
    stripnoise.AddGaussianTailNoise(TestSignal);
}