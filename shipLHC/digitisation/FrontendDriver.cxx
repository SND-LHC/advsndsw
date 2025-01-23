#include "FrontendDriver.h"
#include "StripNoise.h"
#include "InducedCharge.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;


FrontendDriver::FrontendDriver() {}

AdvSignal FrontendDriver::FEDResponse(AdvSignal Signal)
{
    AdvSignal ADCResponse = ADCConversion(Signal); 

    AdvSignal FEDResponseSignal ; 
    std::vector<AdvSignal> temp_FEDResponseSignal;

    StripNoise stripnoise{}; 
    InducedCharge inducedcharge{}; 
    
    if (stripsensor::frontend::ZSModeOption && stripsensor::frontend::NoiseOption)
    {
        FEDResponseSignal = stripnoise.AddGaussianTailNoise(Signal); 
        temp_FEDResponseSignal.push_back(FEDResponseSignal);
        FEDResponseSignal = inducedcharge.Combine(temp_FEDResponseSignal); 
        FEDResponseSignal = ZeroSuppressionAlgorithms(FEDResponseSignal);
    }
    if (!stripsensor::frontend::ZSModeOption)
    {
        FEDResponseSignal = stripnoise.AddGaussianNoise(Signal);
        FEDResponseSignal = stripnoise.AddCMNoise(FEDResponseSignal); 
        temp_FEDResponseSignal.push_back(FEDResponseSignal);
        FEDResponseSignal = inducedcharge.Combine(temp_FEDResponseSignal); 
    }
    return FEDResponseSignal;
}

AdvSignal FrontendDriver::ADCConversion(AdvSignal ResponseSignal)
{
        std::vector<Double_t> NumberofElectrons = ResponseSignal.getIntegratedSignal();
        for (int i = 0; i < NumberofElectrons.size(); i++)
        {
            ADCcount.push_back(stripsensor::frontend::StripNoise + std::ceil(NumberofElectrons[i]/stripsensor::frontend::ElectronperADC));
        }

        AdvSignal ADCResponse(ResponseSignal.getStrips(), ADCcount);

        return ADCResponse; 
}

AdvSignal FrontendDriver::SaturateRange(AdvSignal Signal)
{
    std::vector<Double_t> Charge = Signal.getIntegratedSignal();
    for (int i = 0; i < Charge.size(); i++)
    {
        if (stripsensor::frontend::ZSModeOption)
        {
            if (Charge[i] == 1023)
            {
                Charge[i] = 255; 
            }
            if (Charge[i] > 253)
            {
                Charge[i] = 254; 
            }
        } else {
            if (Charge[i] > 1023)
            {
                Charge[i] = 1023; 
            }
        }
    }
    AdvSignal SaturatedSignal(Signal.getStrips(), Charge);
    return SaturatedSignal; 
}

AdvSignal FrontendDriver::ZeroSuppressionAlgorithms(AdvSignal Signal)
{
    std::vector<Double_t> Amplitude = Signal.getIntegratedSignal(); 
    std::vector<Int_t> Strips = Signal.getStrips();

    std::vector<Int_t> ClusterStrips; 
    std::vector<Double_t> ClusterAmplitudes; 

    std::vector<Int_t> temp_ClusterStrips; 
    std::vector<Double_t> temp_ClusterAmplitudes;

    Int_t mode = stripsensor::frontend::ZeroSuppressionMode; 

    switch (mode) {
        case 1:
            for (int i = 0; i < Amplitude.size(); i++)
            {
                if (Amplitude[i] > stripsensor::frontend::ZeroSuppressionMode1T)
                {
                    ClusterStrips.push_back(Strips[i]);
                    ClusterAmplitudes.push_back(Amplitude[i]);
                }
            } 
            break; 
        case 2:
            for (int i = 0; i < Amplitude.size(); i++)
            {
                if (Amplitude[i] > 2*stripsensor::frontend::StripNoise)
                {
                    temp_ClusterStrips.push_back(Strips[i]); 
                    temp_ClusterAmplitudes.push_back(Amplitude[i]); 
                }
            } 
            for (int j = 0; j < temp_ClusterStrips.size(); j++)
            {
                if((std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+1)==temp_ClusterStrips.end()) && (std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]-1)==temp_ClusterStrips.end())) 
                {
                    if (temp_ClusterAmplitudes[j] > 5*stripsensor::frontend::StripNoise)
                    {
                        ClusterStrips.push_back(temp_ClusterStrips[j]); ClusterAmplitudes.push_back(temp_ClusterAmplitudes[j]);
                    }
                } else {
                    ClusterStrips.push_back(temp_ClusterStrips[j]); ClusterAmplitudes.push_back(temp_ClusterAmplitudes[j]);
                }
            }
            break; 
        case 3:
            for (int i = 0; i < Amplitude.size(); i++)
            {
                if (Amplitude[i] > 3*stripsensor::frontend::StripNoise)
                {
                    temp_ClusterStrips.push_back(Strips[i]); 
                    temp_ClusterAmplitudes.push_back(Amplitude[i]); 
                }
            } 

            for (int j = 0; j < temp_ClusterStrips.size(); j++)
            {
                if(std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+1)!=temp_ClusterStrips.end()) 
                {
                    auto it = std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+1 );
                    if ((temp_ClusterAmplitudes[it - temp_ClusterStrips.begin()])>0)
                    {
                        ClusterStrips.push_back(temp_ClusterStrips[it - temp_ClusterStrips.begin()]);
                        ClusterAmplitudes.push_back(temp_ClusterAmplitudes[it - temp_ClusterStrips.begin()]);
                    }

                } 
                if(std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]-1)!=temp_ClusterStrips.end()) 
                {
                    auto itlower = std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]-1);
                    if ((temp_ClusterAmplitudes[itlower - temp_ClusterStrips.begin()])>0)
                    {
                        ClusterStrips.push_back(temp_ClusterStrips[itlower - temp_ClusterStrips.begin()]);
                        ClusterAmplitudes.push_back(temp_ClusterAmplitudes[itlower - temp_ClusterStrips.begin()]);
                    }

                } 

            }
            break; 
        case 4:
            Double_t SumCharge = 0;
            Double_t SumNoise = 0;

            for (int i = 0; i < Amplitude.size(); i++)
            {
                if (Amplitude[i] > 3*stripsensor::frontend::StripNoise)
                {
                    temp_ClusterStrips.push_back(Strips[i]); 
                    temp_ClusterAmplitudes.push_back(Amplitude[i]); 
                }
            } 
            Int_t neighbhour = 0 ; 
            Int_t neighbhourlower = 0; 
            Int_t j = 0; 
            std::vector<Int_t> Clusters ; 
            while (j < temp_ClusterStrips.size())
            {
                neighbhour = 1; 
                neighbhourlower = 1;
                while((std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+neighbhour)!=temp_ClusterStrips.end())) 
                {
                    if((std::find(ClusterStrips.begin(), ClusterStrips.end(), temp_ClusterStrips[j])==ClusterStrips.end()))
                    {
                        ClusterStrips.push_back(temp_ClusterStrips[j]);
                        ClusterAmplitudes.push_back(temp_ClusterAmplitudes[j]);   
                    }
                    auto it = std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+neighbhour );
                    if ((temp_ClusterAmplitudes[it - temp_ClusterStrips.begin()])>2*stripsensor::frontend::StripNoise)
                    {
                        ClusterStrips.push_back(temp_ClusterStrips[it - temp_ClusterStrips.begin()]);
                        ClusterAmplitudes.push_back(temp_ClusterAmplitudes[it - temp_ClusterStrips.begin()]);
                        neighbhour += 1; 
                    }
                    
                } 
                if (neighbhour > 1)
                {
                    Clusters.push_back(j); 
                    Clusters.push_back(j+neighbhour);
                }
                j = j + neighbhour; 
            }
  
            for (int l = 0; l < Clusters.size(); l+=2)
            {
                for (int m = Clusters[l]; m < Clusters[l+1]; m++)
                {
                    SumCharge += temp_ClusterAmplitudes[m]; 
                    SumNoise += pow(stripsensor::frontend::StripNoise, 2);
                } 
                if (SumCharge < 5*SumNoise)
                {
                     
                    for (int m = Clusters[l]; m < Clusters[l+1]; m++)
                    {
                        ClusterStrips.erase(find(ClusterStrips.begin(), ClusterStrips.end(), temp_ClusterStrips[m]));
                        ClusterAmplitudes.erase(find(ClusterAmplitudes.begin(), ClusterAmplitudes.end(), temp_ClusterAmplitudes[m]));
                    } 
                }
                SumCharge = 0 ; 
                SumNoise = 0 ; 
            }
            break; 
    }
    AdvSignal ClusterSignal(ClusterStrips, ClusterAmplitudes); 
    return ClusterSignal; 
}