#include "sndEventHeaderBaseCut.h"

#include "SNDLHCEventHeader.h"
#include "TChain.h"

#include <stdexcept>

namespace sndAnalysis {

SNDLHCEventHeader* EventHeaderBaseCut::header = 0;
TChain* EventHeaderBaseCut::tree = 0;

EventHeaderBaseCut::EventHeaderBaseCut(TChain* ch)
{
    if (header == 0) {
        header = new SNDLHCEventHeader();
        ch->SetBranchAddress("EventHeader", &header);
        ch->GetEntry(0);
        if (header->GetEventTime() == -1) {
            ch->SetBranchAddress("EventHeader.", &header);
            ch->GetEntry(0);
            if (header->GetEventTime() == -1)
                throw std::runtime_error("Invalid event header");
        }
        tree = ch;
    }
}
}   // namespace sndAnalysis
