#pragma once

#include "sndSciFiBaseCut.h"

#include "TChain.h"
#include "sndScifiHit.h"

namespace sndAnalysis {
  class minSciFiConsecutivePlanes : public sciFiBaseCut {
  public :
    minSciFiConsecutivePlanes(TChain * tree);
    ~minSciFiConsecutivePlanes(){;}

    bool passCut();

  };
};
