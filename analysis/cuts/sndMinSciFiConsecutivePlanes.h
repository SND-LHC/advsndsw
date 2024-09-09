#pragma once

#include "TChain.h"
#include "sndSciFiBaseCut.h"
#include "sndScifiHit.h"

namespace sndAnalysis {
class minSciFiConsecutivePlanes : public sciFiBaseCut
{
  public:
    minSciFiConsecutivePlanes(TChain* tree);
    ~minSciFiConsecutivePlanes() { ; }

    bool passCut();
};
};   // namespace sndAnalysis
