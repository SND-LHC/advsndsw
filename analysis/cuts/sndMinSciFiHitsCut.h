#pragma once

#include "sndSciFiBaseCut.h"

#include "TChain.h"
#include "sndScifiHit.h"

namespace sndAnalysis {
  class minSciFiHits : public sciFiBaseCut {
  private :
    int hitThreshold;
  public :
    minSciFiHits(int threshold, TChain * tree);
    ~minSciFiHits(){;}

    bool passCut();

  };
};
