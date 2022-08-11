#include "SNDLHCEventHeader.h"
#include "SNDLHCEventHeaderConst.h"
#include "FairRootManager.h"

#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// Helper function to get the left shift corresponding to a mask
constexpr int maskToShift(uint64_t mask) {
    int c = 64; // c will be the number of zero bits on the right
    mask &= -mask;
    if (mask) c--;
    if (mask & 0x00000000FFFFFFFF) c -= 32;
    if (mask & 0x0000FFFF0000FFFF) c -= 16;
    if (mask & 0x00FF00FF00FF00FF) c -= 8;
    if (mask & 0x0F0F0F0F0F0F0F0F) c -= 4;
    if (mask & 0x3333333333333333) c -= 2;
    if (mask & 0x5555555555555555) c -= 1;

    return c;
}

// -----   Default constructor   -------------------------------------------
SNDLHCEventHeader::SNDLHCEventHeader()
  : FairEventHeader()
  , fFlags(0)
  , fFillNumber(0)
  , fAccMode(0)
  , fBeamMode(0)
{}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
SNDLHCEventHeader::SNDLHCEventHeader(Int_t runN, uint64_t evtNumber, int64_t timestamp, uint64_t flags)
{
   SetRunId(runN);
   SetEventTime(timestamp);
   SetMCEntryNumber(evtNumber);
   SetFlags(flags);
}

// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
SNDLHCEventHeader::~SNDLHCEventHeader() { }
// -------------------------------------------------------------------------

//-----   Setters   --------------------------------------------------------
void SNDLHCEventHeader::SetFlags(uint64_t flags)
{
    fFlags = flags;
    fFillNumber = (fFlags & FILL_NUMBER_MASK);
    fAccMode    = (fFlags & ACCELERATOR_MODE_MASK) >> maskToShift(ACCELERATOR_MODE_MASK);
    fBeamMode   = (fFlags & BEAM_MODE_MASK) >> maskToShift(BEAM_MODE_MASK);
}
//-----   Getters   ----------------------------------------------------
string SNDLHCEventHeader::GetTimeAsString()
{
   time_t time = fUTCtimestamp;
   struct tm *GMTtime;
   GMTtime = gmtime(&time);

   return asctime(GMTtime);
}

map<string, bool> SNDLHCEventHeader::GetFastNoiseFilters()
{
   map<string, bool> FastNoiseFilters{};

   auto base = (fFlags & FAST_FILTER_MASK) >> maskToShift(FAST_FILTER_MASK);
   FastNoiseFilters["SciFi"]       = base & FAST_FILTER_SCIFI;
   FastNoiseFilters["SciFi_Total"] = base & FAST_FILTER_SCIFI_TOTAL;
   FastNoiseFilters["US"]          = base & FAST_FILTER_US;
   FastNoiseFilters["US_Total"]    = base & FAST_FILTER_US_TOTAL;
   FastNoiseFilters["DS"]          = base & FAST_FILTER_DS;
   FastNoiseFilters["DS_Total"]    = base & FAST_FILTER_DS_TOTAL;
   FastNoiseFilters["Veto_Total"]  = base & FAST_FILTER_VETO_TOTAL;

   return FastNoiseFilters;
}

map<string, bool> SNDLHCEventHeader::GetAdvNoiseFilters()
{
   map<string, bool> AdvNoiseFilters{};

   auto base = (fFlags & ADVANCED_FILTER_MASK) >> maskToShift(ADVANCED_FILTER_MASK);
   AdvNoiseFilters["SciFi_Planes"]  = base & ADVANCED_FILTER_SCIFI_PLANES;
   AdvNoiseFilters["SciFi_Hits"]    = base & ADVANCED_FILTER_SCIFI_HITS;
   AdvNoiseFilters["US_Planes"]     = base & ADVANCED_FILTER_US_PLANES;
   AdvNoiseFilters["US_Hits"]       = base & ADVANCED_FILTER_US_HITS;
   AdvNoiseFilters["DSH_Planes"]    = base & ADVANCED_FILTER_DSH_PLANES;
   AdvNoiseFilters["DSH_Hits"]      = base & ADVANCED_FILTER_DSH_HITS;
   AdvNoiseFilters["DSV_Planes"]    = base & ADVANCED_FILTER_DSV_PLANES;
   AdvNoiseFilters["DSV_Hits"]      = base & ADVANCED_FILTER_DSV_HITS;
   AdvNoiseFilters["DS_Planes"]     = base & ADVANCED_FILTER_DS_PLANES;
   AdvNoiseFilters["VETO_Planes"]   = base & ADVANCED_FILTER_VETO_PLANES;
   AdvNoiseFilters["VETO_Hits"]     = base & ADVANCED_FILTER_VETO_HITS;
   AdvNoiseFilters["GLOBAL_Planes"] = base & ADVANCED_FILTER_GLOBAL_PLANES;

   return AdvNoiseFilters;
}

vector<string> SNDLHCEventHeader::GetPassedFastNFCriteria()
{
  map<string, bool> FastNoiseFilters = GetFastNoiseFilters();
  vector<string> passed;
  for ( auto it : FastNoiseFilters )
      if ( it.second == true ) passed.push_back(it.first);

  return passed;
}

vector<string> SNDLHCEventHeader::GetPassedAdvNFCriteria()
{
  map<string, bool> AdvNoiseFilters = GetAdvNoiseFilters();
  vector<string> passed;
  for ( auto it : AdvNoiseFilters )
      if ( it.second == true ) passed.push_back(it.first);

  return passed;
}

// -----   Override for Fair's ioman->Register   ---------------------------
void SNDLHCEventHeader::Register(bool Persistence)
{
    FairRootManager::Instance()->Register("EventHeader.",
                                          "sndEventHeader",
                                           this, Persistence);
}

// -----   Public method Print   -------------------------------------------
void SNDLHCEventHeader::Print(const Option_t* opt) const
{

  cout << "-I- SNDLHCEventHeader: run number " << fRunId
       << " event number "    << fMCEntryNo 
       << " timestamp "       << fEventTime << endl
       << " LHC fill number " << fFillNumber
       << " LHC beam mode "   << fBeamMode << endl;

}
// -------------------------------------------------------------------------

ClassImp(SNDLHCEventHeader)

