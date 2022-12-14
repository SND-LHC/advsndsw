#pragma once

#include <vector>

#include "sndBaseCut.h"

#include "TChain.h"
#include "TClonesArray.h"

namespace sndAnalysis {
  
  class MuFilterBaseCut : public baseCut {

  protected :
    static TClonesArray * muFilterDigiHitCollection;

    MuFilterBaseCut(TChain * ch);
    ~MuFilterBaseCut(){;}
  };
}
