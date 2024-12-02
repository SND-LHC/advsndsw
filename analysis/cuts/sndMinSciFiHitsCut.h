#pragma once

#include "TChain.h"
#include "sndSciFiBaseCut.h"
#include "sndScifiHit.h"

namespace sndAnalysis {
class minSciFiHits : public sciFiBaseCut
{
  private:
    int hitThreshold;

  public:
    minSciFiHits(int threshold, TChain* tree);
    ~minSciFiHits() { ; }

    bool passCut();
};
};   // namespace sndAnalysis
