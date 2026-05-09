#include "AdvHit.h"
#include "AdvPoint.h"
#include "digitisation/AdvSignal.h"
#include "digitisation/AdvDigitisation.h"
#include "FairLogger.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TStopwatch.h"

#include <TDatabasePDG.h>
#include <iomanip>
#include <typeinfo>
#include <iostream>
#include <string>
#include <map>
using std::cout;
using std::endl;

// -----   Default constructor   -------------------------------------------
AdvHit::AdvHit()
  : TObject()
{
    flag = true;
    signal = -999;
    time = -999;
}
// -----   Standard constructor   ------------------------------------------
AdvHit::AdvHit(Int_t detID)
  :TObject(),
   fDetectorID(detID)
{
    flag = true;
    signal = -999;
    time = -999;
}

// -----   constructor from AdvPoint   ------------------------------------------
AdvHit::AdvHit(Int_t detID, const std::vector<AdvPoint*>& V)
{
    fDetectorID = detID;
    AdvDigitisation advdigi{};
    fDigitisedHit = advdigi.digirunoutput(detID, V);
    flag = true;
    LOG(DEBUG) << "signal created";
}

// -----   Public method Print   -------------------------------------------
void AdvHit::Print() const
{
    LOG(INFO) << " AdvHit: " << " in detector " << fDetectorID
              << "signal " << signal << " time " << time;
}
