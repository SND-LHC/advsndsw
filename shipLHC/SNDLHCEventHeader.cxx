#include "SNDLHCEventHeader.h"
#include "SNDLHCEventHeaderConst.h"
#include "FairRootManager.h"

#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

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
   fFlags = flags;
   SetFillNumber(flags);
   SetBeamMode(flags);
   SetAccMode(flags);
}

// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
SNDLHCEventHeader::~SNDLHCEventHeader() { }
// -------------------------------------------------------------------------

//-----   Setters   --------------------------------------------------------
void SNDLHCEventHeader::SetFlags(uint64_t flags)
{
    fFlags = flags;
    SetFillNumber(flags);
    SetBeamMode(flags);
    SetAccMode(flags);
}

//-----   Getters   ----------------------------------------------------
string SNDLHCEventHeader::GetTimeAsString()
{
   time_t time = fEventTime;
   struct tm *GMTtime;
   GMTtime = gmtime(&time);

   return asctime(GMTtime);
}

map<string, int> SNDLHCEventHeader::GetFastNoiseFilters()
{
   map<string, int> FastNoiseFilters{};
   vector<string> fastNoiseFilters = { "SciFi", "SciFi_Total",
                                        "US", "US_Total",
                                        "DS", "DS_Total",
                                        "Veto_Total" };
   // reading flag bits 26 - 32
   int i = 25;
   for ( auto item : fastNoiseFilters )
   {     
      FastNoiseFilters[item] = ( (fFlags >> i) & 1 );
      i++;
   }
   //for test
   //for(auto it: FastNoiseFilters) cout<<" "<< it.first<< " "<<it.second<<endl;

   return FastNoiseFilters;
}

map<string, int> SNDLHCEventHeader::GetAdvNoiseFilters()
{
   map<string, int> AdvNoiseFilters{};
   vector<string> advNoiseFilters  = { "SciFi_Planes", "SciFi_Hits",
                                        "US_Planes", "US_Hits",
                                        "DSH_Planes", "DSH_Hits",
                                        "DSV_Planes", "DSV_Hits", "DS_Planes",
                                        "Veto_Planes", "Veto_Hits",
                                        "Global_Planes"};
   // reading flag bits 33 - 44
   int i = 32;
   for ( auto item : advNoiseFilters )
   {     
      AdvNoiseFilters[item] = ( (fFlags >> i) & 1 );
      i++;
   }
   //for test
   //for(auto it: AdvNoiseFilters) cout<<" "<< it.first<< " "<<it.second<<endl;

   return AdvNoiseFilters;
}

vector<string> SNDLHCEventHeader::GetPassedFastNFCriteria()
{
  map<string, int> FastNoiseFilters = GetFastNoiseFilters();
  vector<string> passed;
  for ( auto it : FastNoiseFilters )
      if ( it.second >0 ) passed.push_back(it.first);

  return passed;
}

vector<string> SNDLHCEventHeader::GetPassedAdvNFCriteria()
{
  map<string, int> AdvNoiseFilters = GetAdvNoiseFilters();
  vector<string> passed;
  for ( auto it : AdvNoiseFilters )
      if ( it.second >0 ) passed.push_back(it.first);

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
       << " UTC timestamp "   << fEventTime << endl
       << " LHC fill number " << fFillNumber
       << " LHC beam mode "   << fBeamMode << endl;

}
// -------------------------------------------------------------------------

ClassImp(SNDLHCEventHeader)

