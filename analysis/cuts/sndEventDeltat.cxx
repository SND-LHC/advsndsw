#include "sndEventDeltat.h"

#include <iostream>

namespace sndAnalysis {

  eventDeltatCut::eventDeltatCut(int delta_event, int delta_timestamp, TChain * ch) : EventHeaderBaseCut(ch) {
    delta_e = delta_event;
    delta_t = delta_timestamp;

    cutName = std::to_string(delta_e)+" event more than "+std::to_string(delta_t)+" clock cycles away";
  }

  bool eventDeltatCut::passCut(){
    unsigned long int current_entry = tree->GetReadEntry();
    long int current_time = header->GetEventTime();

    bool passes = true;
    tree->GetEntry(current_entry + delta_e);

    int sign = (delta_e > 0) - (delta_e < 0);

    if (-sign*(current_time - header->GetEventTime()) <= delta_t) passes = false;

    // Get current entry back
    tree->GetEntry(current_entry);
    return passes;
  }
}
