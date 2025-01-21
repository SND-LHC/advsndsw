#include "InducedCharge.h"
#include "SurfaceSignal.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TRandom.h"
#include "SiSensor.h"
#include "AdvSignal.h"
#include "SiDigiParameters.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <filesystem>
#include <random>

using namespace std; 

InducedCharge::InducedCharge() {}

std::vector<AdvSignal> InducedCharge::IntegrateCharge(std::vector<SurfaceSignal> DiffusionSignal)
{   
    std::vector<AdvSignal> ResponseSignal; 
    for (int k = 0; k < DiffusionSignal.size(); k++)
    {
        std::vector<Double_t> amplitude = DiffusionSignal[k].getAmplitude(); 
        std::vector<TVector3> surfacepos = DiffusionSignal[k].getSurfacePos();
        std::vector<Double_t> diffusionarea = DiffusionSignal[k].getDiffusionArea();

        std::vector<Int_t> AffectedStrips;
        std::vector<Int_t> temp_AffectedStrips; 
        std::vector<Double_t> ChargeDeposited; 
        std::vector<Double_t> TotalChargeDeposited;
        std::vector<std::vector<Double_t>> PulseResponse; 

        Int_t N ;
        Double_t x_start; 
        Double_t x_end; 
        Double_t z_start; 
        Double_t z_end;  
        Double_t integratedcharge; 

        Double_t e = 1.6e-19; 

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
            TotalChargeDeposited[k] = (TotalChargeDeposited[k] * rescale_ratio)*e ; 
        }

        PulseResponse = stripsensor::peakmode ? GetPulseShape(stripsensor::APVpeakpulse, TotalChargeDeposited): GetPulseShape(stripsensor::APVdecopulse, TotalChargeDeposited) ;
        AdvSignal PulseSignal(UniqueAffectedStrips, TotalChargeDeposited, PulseResponse); 
        ResponseSignal.push_back(PulseSignal); 
    }
    return ResponseSignal;
}

std::vector<Int_t> InducedCharge::GetStrips(TVector3 point, Double_t area)
{
    std::vector<Int_t> affectedstrips;

    int fromstrip = floor(((point.X()-(NSigma*area)) / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2)); // check calculation 
    fromstrip = max(0, fromstrip);
    fromstrip = min(advsnd::strips - 1, fromstrip); 

    int tostrip = floor(((point.X()+(NSigma*area)) / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2));
    tostrip = max(0, tostrip);
    tostrip = min(advsnd::strips - 1, tostrip);

    Int_t N; 
    N = tostrip - fromstrip; 

    for (int i = 0 ; i <= N ; i++)
    {
        affectedstrips.push_back(fromstrip + i);
    }

    return affectedstrips; 
}

std::vector<std::vector<Double_t>> InducedCharge::GetPulseShape(std::string PulseFileName, std::vector<Double_t> ChargeDeposited)
{
    std::vector<double> PulseValues;

    std::ifstream inputFile(PulseFileName);

    if (!inputFile.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
    }
    std::string line;
    std::string res_find = "resolution:";
    float res;
    std::string s;

    while (getline(inputFile, line)) {
        if ((!line.empty()) && (line.substr(0, 1) != "#")) {
            std::stringstream ss(line);
            if (line.find(res_find) != std::string::npos) {
                res = stof(line.substr(line.find(res_find) + res_find.length()));   // implement check for
                                                                                    // resolution
            } else {
                std::string value;
                while (getline(ss, value, ' ')) {
                    PulseValues.push_back(stod(value));
                }
            }
        }
    }
    const auto max_value = max_element(PulseValues.begin(), PulseValues.end());
    if (abs(*max_value - 1) > numeric_limits<double>::epsilon()) {
        throw invalid_argument("Maximum value of pulse shape not 1.");
    }

    unsigned int pulset0Idx = std::distance(PulseValues.begin(), max_value);

    std::vector<std::vector<Double_t>> PulseResponse; 
    std::vector<Double_t> temp_response ;
    Double_t response_value; 

    for(int i = 0; i < ChargeDeposited.size(); i++)
    {
         
        Double_t amplitude_max = ChargeDeposited[i]; 

        std::default_random_engine generator;
        std::normal_distribution<double> dist(stripsensor::noise_mean, stripsensor::noise_std_dev);

        for(int j = 0; j < PulseValues.size(); j++)
        {  
            response_value = ((stripsensor::baseline + (PulseValues[j] * amplitude_max * stripsensor::amplificaton_factor)) > stripsensor::rail ? stripsensor::rail : (stripsensor::baseline + (PulseValues[j] * amplitude_max * stripsensor::amplificaton_factor)));

            temp_response.push_back(response_value + dist(generator));
        }

        PulseResponse.push_back(temp_response);
    } 
    return PulseResponse; 
    // get the vector of beginning to max of the pulse
    // time response not included!
}