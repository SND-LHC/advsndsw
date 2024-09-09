#pragma once

#include "TChain.h"
#include "sndMuFilterBaseCut.h"

namespace sndAnalysis {

class vetoCut : public MuFilterBaseCut
{
  public:
    vetoCut(TChain* ch);
    ~vetoCut() { ; }
    bool passCut();
};
}   // namespace sndAnalysis
