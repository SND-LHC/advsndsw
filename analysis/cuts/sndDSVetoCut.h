#pragma once

#include "TChain.h"
#include "sndMuFilterBaseCut.h"

namespace sndAnalysis {

class DSVetoCut : public MuFilterBaseCut
{
  public:
    DSVetoCut(TChain* ch);
    ~DSVetoCut() { ; }
    bool passCut();
};
}   // namespace sndAnalysis
