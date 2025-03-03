#include "SndlhcHit.h"

// -----   Default constructor   -------------------------------------------
SndlhcHit::SndlhcHit()
    : TObject()
    , fDetectorID(-1)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
SndlhcHit::SndlhcHit(Int_t detID, Int_t nP, Int_t nS)
    : TObject()
    , fDetectorID(detID)
    , nSiPMs(nP)
    , nSides(nS)
{
    for (unsigned int j = 0; j < 16; ++j) {
        signals[j] = -1;
        times[j] = -1;
        fDaqID[j] = -1;
    }
}

Float_t SndlhcHit::GetSignal(Int_t nChannel) { return signals[nChannel]; }
Float_t SndlhcHit::GetTime(Int_t nChannel) { return times[nChannel]; }
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
SndlhcHit::~SndlhcHit() {}
// -------------------------------------------------------------------------

ClassImp(SndlhcHit)
