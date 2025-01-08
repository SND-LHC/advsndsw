#include <iostream>
#include <algorithm>

#include "Clustering.h"
#include "SiDigiParameters.h"
#include "AdvSignal.h"

Clustering::Clustering() {}

void Clustering::FindClusters(std::vector<AdvSignal> ResponseSignal)
{

    for (int k = 0; k < ResponseSignal.size(); k++)
    {
        std::vector<Double_t> ADCValues = ResponseSignal[k].getIntegratedSignal();
        Double_t MaxADC = *std::max_element(ADCValues.begin(), ADCValues.end());
        if (MaxADC > 3*stripsensor::frontend::StripNoise) 
            {
                std::cout << MaxADC << endl; 
            }
    }
}