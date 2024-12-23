#pragma once

#include "TChain.h"
#include "sndEventHeaderBaseCut.h"

namespace sndAnalysis {

class eventDeltatCut : public EventHeaderBaseCut
{
  private:
    int delta_t;
    int delta_e;

  public:
    eventDeltatCut(int delta_event, int delta_timestamp, TChain* ch);
    ~eventDeltatCut() { ; }
    bool passCut();
};

}   // namespace sndAnalysis
