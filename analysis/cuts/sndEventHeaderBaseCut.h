#pragma once

#include "SNDLHCEventHeader.h"
#include "TChain.h"
#include "sndBaseCut.h"

#include <vector>

namespace sndAnalysis {

class EventHeaderBaseCut : public baseCut
{

  protected:
    static SNDLHCEventHeader* header;
    static TChain* tree;

    EventHeaderBaseCut(TChain* ch);
    ~EventHeaderBaseCut() { ; }
};
}   // namespace sndAnalysis
