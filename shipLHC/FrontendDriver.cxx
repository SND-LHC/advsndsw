#include "FrontendDriver.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;


FrontendDriver::FrontendDriver() {}

// std::vector<AdvSignal> FrontendDriver::ADCConversion(std::vector<AdvSignal> ResponseSignal)
// {
//     for (int k = 0; k < ResponseSignal.size(); k++)
//     {
//         std::vector<Double_t> NumberofElectrons = ResponseSignal[k].getIntegratedSignal();
//         for (int i = 0; i < NumberofElectrons.size(); i++)
//         {
//             ADCcount.push_back(stripsensor::frontend::StripNoise + std::ceil(NumberofElectrons[i]/stripsensor::frontend::ElectronperADC));
//         }

//         AdvSignal ADCResponse(ResponseSignal[k].getStrips(), ADCcount);
//         FEDResponse.push_back(ADCResponse);
//         //ZeroSuppressionAlgorithms(stripsensor::frontend::ZeroSuppressionMode);
//     }
//     return FEDResponse; 
// }

// void FrontendDriver::ZeroSuppressionAlgorithms(AdvSignal Signal)
// {
//     std::vector<Double_t> Amplitude = Signal.getIntegratedSignal(); 
//     std::vector<Int_t> Strips = Signal.getStrips();

//     std::vector<Int_t> ClusterStrips; 
//     std::vector<Double_t> ClusterAmplitudes; 

//     std::vector<Int_t> temp_ClusterStrips; 
//     std::vector<Double_t> temp_ClusterAmplitudes;

//     Int_t mode = stripsensor::frontend::ZeroSuppressionMode; 

//     switch (mode) {
//         case 1:
//             for (int i = 0; i < Amplitude.size(); i++)
//             {
//                 if (Amplitude[i] > stripsensor::frontend::ZeroSuppressionMode1T)
//                 {
//                     ClusterStrips.push_back(Strips[i]);
//                     ClusterAmplitudes.push_back(Amplitude[i]);
//                     cout << Amplitude[i] << endl; 
//                 }
//             } 
//             break; 
//         case 2:
//             for (int i = 0; i < Amplitude.size(); i++)
//             {
//                 if (Amplitude[i] > 2*stripsensor::frontend::StripNoise)
//                 {
//                     temp_ClusterStrips.push_back(Strips[i]); 
//                     temp_ClusterAmplitudes.push_back(Amplitude[i]); 
//                 }
//             } 
//             if (temp_ClusterStrips.size() == 1) 
//             {
//                 if (temp_ClusterStrips[0] > 5*stripsensor::frontend::StripNoise){ClusterStrips = temp_ClusterStrips; ClusterAmplitudes = temp_ClusterAmplitudes;}
//                 else {ClusterStrips = temp_ClusterStrips; ClusterAmplitudes.push_back(0);}
//             }
//             else {ClusterStrips = temp_ClusterStrips; ClusterAmplitudes = temp_ClusterAmplitudes;}
//             for (int i = 0; i < ClusterAmplitudes.size(); i++)
//             {
//                 cout << ClusterStrips[i] << "\t" << ClusterAmplitudes[i] << endl; 
//             } 
//             break; 
//         case 3:
//             for (int i = 0; i < Amplitude.size(); i++)
//             {
//                 if (Amplitude[i] > 3*stripsensor::frontend::StripNoise)
//                 {
//                     temp_ClusterStrips.push_back(Strips[i]); 
//                     temp_ClusterAmplitudes.push_back(Amplitude[i]); 
//                 }
//             } 
//             if (Strips.size() > temp_ClusterStrips.size() && temp_ClusterStrips.size()!=0 ) 
//             {
//                 if(std::find(Strips.begin(), Strips.end(), temp_ClusterStrips[0]+1)!=Strips.end()){
//                 auto it = std::find(Strips.begin(), Strips.end(), temp_ClusterStrips[0]+1);
//                     if ((Amplitude[it - Strips.begin()])>0)
//                     {
//                         temp_ClusterStrips.push_back(Strips[it - Strips.begin()]);
//                         temp_ClusterAmplitudes.push_back(Amplitude[it - Strips.begin()]);
//                     } else {
//                         temp_ClusterAmplitudes.push_back(0);
//                     }
//                 }
//                 if(std::find(Strips.begin(), Strips.end(), temp_ClusterStrips[0]-1)!=Strips.end()){
//                 auto it = std::find(Strips.begin(), Strips.end(), temp_ClusterStrips[0]-1);
//                     if ((Amplitude[it - Strips.begin()])>0)
//                     {
//                         temp_ClusterStrips.emplace(temp_ClusterStrips.begin(), Strips[it - Strips.begin()]);
//                         temp_ClusterAmplitudes.emplace(temp_ClusterAmplitudes.begin(), Amplitude[it - Strips.begin()]);
//                     } else {
//                         temp_ClusterAmplitudes.emplace(temp_ClusterAmplitudes.begin(), 0);
//                     }
//                 }
//             }
//             ClusterStrips = temp_ClusterStrips; 
//             ClusterAmplitudes = temp_ClusterAmplitudes; 
//             for(int k = 0; k < ClusterStrips.size(); k++)
//             {
//                 cout << ClusterStrips[k] << "\t" << ClusterAmplitudes[k] << endl;
//             }
//             break; 
//         case 4:
//             Double_t SumCharge = 0;
//             Double_t SumNoise = 0;
//             for (int i = 0; i < Amplitude.size(); i++)
//             {
//                 if (Amplitude[i] > 3*stripsensor::frontend::StripNoise)
//                 {
//                     temp_ClusterStrips.push_back(Strips[i]); 
//                     temp_ClusterAmplitudes.push_back(Amplitude[i]); 
//                     SumCharge += Amplitude[i];
//                     SumNoise += pow((stripsensor::frontend::StripNoise), 2);
//                 } else if (Amplitude[i] > 2*stripsensor::frontend::StripNoise) {
//                     SumCharge += Amplitude[i];
//                     SumNoise += pow((stripsensor::frontend::StripNoise), 2);
//                     cout << SumCharge << "\t" << SumNoise << endl; 
//                     if (SumCharge > 5*SumNoise)
//                     {
//                         temp_ClusterStrips.push_back(Strips[i]); 
//                         temp_ClusterAmplitudes.push_back(Amplitude[i]);
//                     }
//                 }
//             } 
//             ClusterStrips = temp_ClusterStrips; 
//             ClusterAmplitudes = temp_ClusterAmplitudes;
//             for (int k = 0; k < ClusterStrips.size(); k++)
//             {
//                 cout << ClusterStrips[k] << "\t" << ClusterAmplitudes[k] << endl; 
//             }
//             break; 
//     }
// }

// void TestingAlgorithm()
// {
//     std::vector<Int_t> Strips = {100, 101, 102}; 
//     std::vector<Double_t> IntegratedSignal = {200, 600, 50};
    
//     AdvSignal TestSignal(Strips, IntegratedSignal); 

//     FrontendDriver frontenddriver;
//     frontenddriver.ZeroSuppressionAlgorithms(TestSignal);

// }