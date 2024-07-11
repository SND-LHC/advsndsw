#ifndef SNDFAIRTASKS_DIGITASKSND_H_
#define SNDFAIRTASKS_DIGITASKSND_H_

#include "FairMCEventHeader.h"   // for FairMCEventHeader
#include "FairTask.h"            // for FairTask, InitStatus
#include "SNDLHCEventHeader.h"   // for EventHeader
#include "Scifi.h"               // for Scifi detector

#include <Rtypes.h>       // for THashConsistencyHolder, ClassDef
#include <RtypesCore.h>   // for Double_t, Int_t, Option_t
#include <TClonesArray.h>
class TBuffer;
class TClass;
class TClonesArray;
class TMemberInspector;

using namespace std;

class DigiTaskSND : public FairTask
{
  public:
    /** Default constructor **/
    DigiTaskSND();

    /** Destructor **/
    ~DigiTaskSND();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

  private:
    void digitizeMuFilter();
    void digitizeScifi();
    void digitiseAdvTarget();
    void digitiseAdvMuFilter();

    Scifi* scifi;
    map<Int_t, map<Int_t, array<float, 2>>> fibresSiPM;
    map<Int_t, map<Int_t, array<float, 2>>> siPMFibres;

    // Input
    FairMCEventHeader* fMCEventHeader;
    TClonesArray* fMuFilterPointArray;   // MC points
    TClonesArray* fScifiPointArray;
    TClonesArray* AdvTargetPoints;
    TClonesArray* AdvMuFilterPoints;
    // Output
    SNDLHCEventHeader* fEventHeader;
    TClonesArray* fMuFilterDigiHitArray;   // hit class (digitized!)
    TClonesArray* fScifiDigiHitArray;
    TClonesArray* AdvTargetHits;
    TClonesArray* AdvMuFilterHits;
    TClonesArray* fMuFilterHit2MCPointsArray;   // link to MC truth
    TClonesArray* fScifiHit2MCPointsArray;
    TClonesArray* AdvTargetHits2MCPoints;
    TClonesArray* AdvMuFilterHits2MCPoints;
    TClonesArray* fvetoPointArray;
    TClonesArray* fEmulsionPointArray;

    TClonesArray* fMCTrackArray;
    DigiTaskSND(const DigiTaskSND&);
    DigiTaskSND& operator=(const DigiTaskSND&);

    ClassDef(DigiTaskSND, 3);
};

#endif   // SNDFAIRTASKS_DIGITASKSND_H_
