#include "FrontendDriver.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;


FrontedDriver::FrontedDriver() {}

std::vector<AdvSignal> FrontedDriver::ADCConversion(std::vector<AdvSignal> ResponseSignal)
{
    for (int k = 0; k < ResponseSignal.size(); k++)
    {
        std::vector<Double_t> NumberofElectrons = ResponseSignal[k].getIntegratedSignal();
        for (int i = 0; i < NumberofElectrons.size(); i++)
        {
            ADCcount.push_back(stripsensor::frontend::StripNoise + std::ceil(NumberofElectrons[i]/stripsensor::frontend::ElectronperADC));
        }

        AdvSignal ADCResponse(ResponseSignal[k].getStrips(), ADCcount);
        FEDResponse.push_back(ADCResponse);
        //ZeroSuppressionAlgorithms(stripsensor::frontend::ZeroSuppressionMode);
    }
    return FEDResponse; 
}

void FrontedDriver::ZeroSuppressionAlgorithms(Int_t mode, AdvSignal Signal)
{
    std::vector<Double_t> Amplitude = Signal.getIntegratedSignal(); 
    std::vector<Int_t> Strips = Signal.getStrips();

    std::vector<Int_t> ClusterStrips; 
    std::vector<Double_t> ClusterAmplitudes; 

    std::vector<Int_t> temp_ClusterStrips; 
    std::vector<Double_t> temp_ClusterAmplitudes;

    Int_t count = 0; 

    switch (mode) {
        case 1:
            for (int i = 0; i < Amplitude.size(); i++)
            {
                if (Amplitude[i] > stripsensor::frontend::ZeroSuppressionMode1T)
                {
                    ClusterStrips.push_back(Strips[i]);
                    ClusterAmplitudes.push_back(Amplitude[i]);
                    cout << Amplitude[i] << endl; 
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
            if (ClusterStrips.size() == 1) 
            {
                if (temp_ClusterStrips[0] > 5*stripsensor::frontend::StripNoise){ClusterStrips = temp_ClusterStrips; ClusterAmplitudes = temp_ClusterAmplitudes;}
                else {ClusterStrips = temp_ClusterStrips; ClusterAmplitudes[0] = 0;}
            }
            break; 
        case 3:
            cout << "Case 3" << endl; 
            break; 
        case 4:
            cout << "Case 4" << endl; 
            break; 
    }
}