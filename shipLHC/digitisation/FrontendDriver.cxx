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

AdvSignal FrontendDriver::FEDResponse(Int_t detID)
{
    std::vector<AdvSignal> temp_FEDResponseSignal;

    StripNoise stripnoise{}; 
    InducedCharge inducedcharge{}; 

    //AdvSignal noisesignal = stripnoise.AddGaussianNoise(detID); 
    AdvSignal pedestalsignal = stripnoise.AddPedestals(detID);
    //temp_FEDResponseSignal.push_back(noisesignal); 
    temp_FEDResponseSignal.push_back(pedestalsignal); 
    AdvSignal FEDResponseSignal = inducedcharge.Combine(temp_FEDResponseSignal); 
    FEDResponseSignal = ZeroSuppressionAlgorithms(detID, FEDResponseSignal);
    FEDResponseSignal = SaturateRange(FEDResponseSignal);
    return FEDResponseSignal; 
}

AdvSignal FrontendDriver::FEDResponse(Int_t detID, AdvSignal Signal)
{
    AdvSignal ADCResponse = ADCConversion(Signal); 

    AdvSignal FEDResponseSignal ; 
    std::vector<AdvSignal> temp_FEDResponseSignal;

    StripNoise stripnoise{}; 
    InducedCharge inducedcharge{}; 
    
    if (stripsensor::frontend::ZSModeOption)
    {
        //AdvSignal noisesignal = stripnoise.AddGaussianNoise(detID);
        AdvSignal pedestalsignal = stripnoise.AddPedestals(detID);
        FEDResponseSignal = stripnoise.AddGaussianTailNoise(Signal); 
        temp_FEDResponseSignal.push_back(FEDResponseSignal);
        temp_FEDResponseSignal.push_back(pedestalsignal);
        //temp_FEDResponseSignal.push_back(noisesignal);
        FEDResponseSignal = inducedcharge.Combine(temp_FEDResponseSignal); 
        FEDResponseSignal = ZeroSuppressionAlgorithms(detID, FEDResponseSignal);
    }
    if (!stripsensor::frontend::ZSModeOption)
    {
        //FEDResponseSignal = stripnoise.AddGaussianNoise(Signal);
        FEDResponseSignal = stripnoise.AddCMNoise(FEDResponseSignal); 
        temp_FEDResponseSignal.push_back(FEDResponseSignal);
        FEDResponseSignal = inducedcharge.Combine(temp_FEDResponseSignal); 
    }
    FEDResponseSignal = SaturateRange(FEDResponseSignal);
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
            if (Charge[i] > 1022)
            {
                Charge[i] = 255; 
            }
            else if (Charge[i] > 253)
            {
                Charge[i] = 254; 
            }
            if (Charge[i] < 0)
            {
                Charge[i] = 0; 
            }
   
        } else {
            if (Charge[i] > 1023)
            {
                Charge[i] = 1023; 
            }
            if (Charge[i] < 0)
            {
                Charge[i] = 0; 
            }
        }
    }
    AdvSignal SaturatedSignal(Signal.getStrips(), Charge);
    return SaturatedSignal; 
}

AdvSignal FrontendDriver::ZeroSuppressionAlgorithms(Int_t detID, AdvSignal Signal)
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
                Int_t current_detid = ((detID & ~0x3FF) | (Strips[i] & 0x3FF)); 
                if (Amplitude[i] > stripsensor::frontend::ZeroSuppressionMode1T*StripNoise::noise_map[current_detid])
                {
                    ClusterStrips.push_back(Strips[i]);
                    ClusterAmplitudes.push_back(Amplitude[i]);
                }
            } 
            break; 
        case 2:
            for (int i = 0; i < Amplitude.size(); i++)
            {
                Int_t current_detid = ((detID & ~0x3FF) | (Strips[i] & 0x3FF)); 
                if (Amplitude[i] > 2*StripNoise::noise_map[current_detid])
                {
                    temp_ClusterStrips.push_back(Strips[i]); 
                    temp_ClusterAmplitudes.push_back(Amplitude[i]); 
                }
            } 

            for (int j = 0; j < temp_ClusterStrips.size(); j++)
            {
                if((std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+1)==temp_ClusterStrips.end()) && (std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]-1)==temp_ClusterStrips.end())) 
                {
                    Int_t current_detid = ((detID & ~0x3FF) | (temp_ClusterStrips[j] & 0x3FF)); 
                    if (temp_ClusterAmplitudes[j] > 5*StripNoise::noise_map[current_detid])
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
                Int_t current_detid = ((detID & ~0x3FF) | (Strips[i] & 0x3FF));
                if (Amplitude[i] > 3*StripNoise::noise_map[current_detid])
                {
                    temp_ClusterStrips.push_back(Strips[i]); 
                    temp_ClusterAmplitudes.push_back(Amplitude[i]); 
                }
            } 

            ClusterStrips = temp_ClusterStrips; 
            ClusterAmplitudes = temp_ClusterAmplitudes; 

            for (int j = 0; j < temp_ClusterStrips.size(); j++)
            {
                if((std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+1)!=temp_ClusterStrips.end()) && (std::find(ClusterStrips.begin(), ClusterStrips.end(), temp_ClusterStrips[j]+1)==ClusterStrips.end()))
                {
                    auto it = std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]+1 );
                    if ((temp_ClusterAmplitudes[it - temp_ClusterStrips.begin()])>0)
                    {
                        ClusterStrips.push_back(temp_ClusterStrips[it - temp_ClusterStrips.begin()]);
                        ClusterAmplitudes.push_back(temp_ClusterAmplitudes[it - temp_ClusterStrips.begin()]);
                    }

                } 
                if((std::find(temp_ClusterStrips.begin(), temp_ClusterStrips.end(), temp_ClusterStrips[j]-1)!=temp_ClusterStrips.end()) && (std::find(ClusterStrips.begin(), ClusterStrips.end(), temp_ClusterStrips[j]-1)==ClusterStrips.end()))
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
            std::vector<int> mask; 
            std::vector<int> h_thr; 
            std::vector<int> l_thr; 

            for (int i =0; i < Strips.size(); i++)
            {
                Int_t current_detid = ((detID & ~0x3FF) | (Strips[i] & 0x3FF)); 
                h_thr.push_back(stripsensor::frontend::hthr * StripNoise::noise_map[current_detid]); 
                l_thr.push_back(stripsensor::frontend::lthr * StripNoise::noise_map[current_detid]);
                mask.push_back(0); 
            }
            for (int j = 0; j < Strips.size(); j++)
            {
                if (Amplitude[j] >= h_thr[j])
                {
                    mask[j] = 1; 
                }
                else if (Amplitude[j] >= l_thr[j] && Amplitude[j] <h_thr[j])  
                {
                    mask[j] = 2; 
                }
                else 
                {
                    mask[j] = 0; 
                }

            }
            for (int k = 0; k < mask.size(); k++)
            {
                if (mask[k] == 2)
                {
                    if (((mask[k-1] == 2 || mask[k-1] == 1) && (abs(Strips[k] - Strips[k-1]) == 1)) || ((mask[k+1] == 2 || mask[k+1] == 1) && (abs(Strips[k] - Strips[k+1]) == 1)))
                    {
                        mask[k] = 1; 
                    }
                }
                if (mask[k] == 0)
                {
                    if ((mask[k-1]==1 && mask[k+1]==1) && abs((Strips[k] - Strips[k-1]))==1 && abs((Strips[k] - Strips[k+1]))==1)
                    {
                        mask[k] = 1; 
                    }
                    else if (((mask[k+1]==1) && (abs((Strips[k] - Strips[k+1]))==1)) && ((mask[k-1]==2) && (abs((Strips[k] - Strips[k-1]))==1)) && ((mask[k-2]==2) && (abs((Strips[k] - Strips[k-2]))==1)))
                    {
                        mask[k] = 1; 
                    }
                    else if (((mask[k-1]==1) && (abs((Strips[k] - Strips[k-1]))==1)) && ((mask[k+1]==2) && (abs((Strips[k] - Strips[k+1]))==1)) && ((mask[k+2]==2) && (abs((Strips[k] - Strips[k+2]))==2)))
                    {
                        mask[k] = 1; 
                    }
                    else if (((mask[k-2]==2) && (abs((Strips[k] - Strips[k-2]))==2)) && ((mask[k-1]==2) && (abs((Strips[k] - Strips[k-1]))==1)) && ((mask[k+1]==2) && (abs((Strips[k] - Strips[k+1]))==1)) && ((mask[k+2]==2) && (abs((Strips[k] - Strips[k+2]))==2)))
                    {
                        mask[k] = 1; 
                    }

                }
            }
            for (int n = 0; n < mask.size(); n++)
            {
                if (mask[n] == 1)
                {
                    ClusterStrips.push_back(Strips[n]); 
                    ClusterAmplitudes.push_back(Amplitude[n]); 
                }
            }

            break; 
    }
    AdvSignal ClusterSignal(ClusterStrips, ClusterAmplitudes); 
    return ClusterSignal; 
}
