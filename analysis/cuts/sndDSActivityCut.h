#pragma once

#include "TChain.h"
#include "sndMuFilterBaseCut.h"

namespace sndAnalysis {

class DSActivityCut : public MuFilterBaseCut
{
  public:
    DSActivityCut(TChain* ch);
    ~DSActivityCut() { ; }
    bool passCut();
};
}   // namespace sndAnalysis
