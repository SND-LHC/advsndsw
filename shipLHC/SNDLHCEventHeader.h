#ifndef SNDLHCEVENTHEADER_H
#define SNDLHCEVENTHEADER_H 1

#include "FairEventHeader.h"
#include "SNDLHCEventHeaderConst.h"

#include <ctime>
#include <map>
#include <string>
#include <vector>

using namespace std;

class SNDLHCEventHeader : public FairEventHeader
{

  public:

    /** Default constructor **/
    SNDLHCEventHeader();

    /** Constructor with arguments **/
    SNDLHCEventHeader(Int_t runN, uint64_t evtNumber, int64_t timestamp, uint64_t flags);

    /** Destructor **/
    virtual ~SNDLHCEventHeader();

    /** Setters **/
    void SetFlags(uint64_t flags);
    void SetFillNumber(uint64_t flags) { fFillNumber = (flags & FILL_NUMBER_MASK); }
    void SetAccMode(uint64_t flags) { fAccMode = static_cast<int>(flags & ACCELERATOR_MODE_MASK); }
    void SetBeamMode(uint64_t flags) { fBeamMode = static_cast<int>(flags & BEAM_MODE_MASK); }

    /** Getters **/
    string GetTimeAsString(); // GMT time
    uint16_t GetFillNumber() const { return fFillNumber; }
    int GetAccMode() const { return fAccMode; }
    int GetBeamMode() const { return fBeamMode; }
    map<string, int> GetFastNoiseFilters();
    map<string, int> GetAdvNoiseFilters();
    vector<string> GetPassedFastNFCriteria();
    vector<string> GetPassedAdvNFCriteria();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    // Override FairEventHeader's Register
    virtual void Register(Bool_t Persistance = kTRUE);

  private:
    uint64_t fFlags;
    uint16_t fFillNumber;
    int fAccMode; // enum class
    int fBeamMode; // enum class
    
    /** Copy constructor **/
    SNDLHCEventHeader(const SNDLHCEventHeader& eventHeader);
    SNDLHCEventHeader operator=(const SNDLHCEventHeader& eventHeader);

    ClassDef(SNDLHCEventHeader,1)

};

#endif /* SNDLHCEVENTHEADER_H */
