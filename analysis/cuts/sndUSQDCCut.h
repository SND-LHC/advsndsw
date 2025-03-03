#pragma once

#include "TChain.h"
#include "sndMuFilterBaseCut.h"

namespace sndAnalysis {

class USQDCCut : public MuFilterBaseCut
{
  private:
    float qdc_threshold;

  public:
    USQDCCut(float threshold, TChain* ch);
    ~USQDCCut() { ; }
    bool passCut();
};
}   // namespace sndAnalysis
