#pragma once

#include "TChain.h"
#include "sndSciFiBaseCut.h"
#include "sndScifiHit.h"

namespace sndAnalysis {
class sciFiStationCut : public sciFiBaseCut
{
  private:
    float fractionThreshold;
    std::vector<int> stations_to_exclude;

  public:
    sciFiStationCut(float threshold, std::vector<int> excluded_stations, TChain* tree);
    ~sciFiStationCut() { ; }

    bool passCut();
};
};   // namespace sndAnalysis
