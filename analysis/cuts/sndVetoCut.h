#pragma once

#include "sndMuFilterBaseCut.h"

#include "TChain.h"

namespace sndAnalysis {
  
  class vetoCut : public MuFilterBaseCut {
  public :
    vetoCut(TChain * ch);
    ~vetoCut(){;}
    bool passCut();
  };
}
