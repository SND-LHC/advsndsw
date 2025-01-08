#include <iostream>
#include <vector>

#include "AdvSignal.h"


class Clustering 
{
    public:
        Clustering();
        void FindClusters(std::vector<AdvSignal> ResponseSignal);

    private:
        std::vector<Double_t> HitPosition; 
        std::vector<Int_t> ClusterStrips;
};