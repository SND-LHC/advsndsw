#pragma once

#include "sndMuFilterBaseCut.h"

#include "TChain.h"

namespace sndAnalysis {
  
  class DSActivityCut : public MuFilterBaseCut {
  public :
    DSActivityCut(TChain * ch);
    ~DSActivityCut(){;}
    bool passCut();
  };
}
