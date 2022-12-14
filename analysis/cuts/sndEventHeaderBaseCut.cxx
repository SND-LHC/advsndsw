#include "sndEventHeaderBaseCut.h"

#include "SNDLHCEventHeader.h"
#include "TChain.h"

namespace sndAnalysis {

  SNDLHCEventHeader * EventHeaderBaseCut::header = 0;
  TChain * EventHeaderBaseCut::tree = 0;

  EventHeaderBaseCut::EventHeaderBaseCut(TChain * ch){
    if (header == 0){
      header = new SNDLHCEventHeader();
      ch->SetBranchAddress("EventHeader", &header);
      tree = ch;
    }
  }
}
