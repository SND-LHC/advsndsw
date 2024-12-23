#pragma once

#include "TChain.h"
#include "TClonesArray.h"
#include "sndBaseCut.h"

#include <vector>

namespace sndAnalysis {

class MuFilterBaseCut : public baseCut
{

  protected:
    static TClonesArray* muFilterDigiHitCollection;

    MuFilterBaseCut(TChain* ch);
    ~MuFilterBaseCut() { ; }
};
}   // namespace sndAnalysis
