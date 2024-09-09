#ifndef Scifi_H
#define Scifi_H

#include "FairDetector.h"
#include "FairModule.h"   // for FairModule
#include "Rtypes.h"       // for ShipMuonShield::Class, Bool_t, etc
#include "SNDLHCEventHeader.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <string>   // for string

class ScifiPoint;
class FairVolume;
class TClonesArray;

class Scifi : public FairDetector
{
  public:
    Scifi(const char* name, Bool_t Active, const char* Title = "Scifi");
    Scifi();
    virtual ~Scifi();

    /**      Create the detector geometry        */
    void ConstructGeometry();

    /** Get position of single fibre in global coordinate system**/
    void GetPosition(Int_t id, TVector3& vLeft, TVector3& vRight);   // or top and bottom
    /** Transform global position to local position in plane **/
    TVector3 GetLocalPos(Int_t id, TVector3* glob);
    /** mean position of fibre2 associated with SiPM channel **/
    void GetSiPMPosition(Int_t SiPMChan, TVector3& A, TVector3& B);
    Double_t GetCorrectedTime(Int_t fDetectorID, Double_t rawTime, Double_t L);
    Double_t ycross(Double_t a, Double_t R, Double_t x);
    Double_t integralSqrt(Double_t ynorm);
    Double_t fraction(Double_t R, Double_t x, Double_t y);
    Double_t area(Double_t a, Double_t R, Double_t xL, Double_t xR);
    void SiPMmapping();
    std::map<Int_t, std::map<Int_t, std::array<float, 2>>> GetSiPMmap() { return fibresSiPM; }
    std::map<Int_t, std::map<Int_t, std::array<float, 2>>> GetFibresMap() { return siPMFibres; }
    std::map<Int_t, float> GetSiPMPos() { return SiPMPos; }
    virtual void SiPMOverlap();
    void SetConfPar(TString name, Float_t value) { conf_floats[name] = value; }
    void SetConfPar(TString name, Int_t value) { conf_ints[name] = value; }
    void SetConfPar(TString name, TString value) { conf_strings[name] = value; }
    Float_t GetConfParF(TString name) { return conf_floats[name]; }
    Int_t GetConfParI(TString name) { return conf_ints[name]; }
    TString GetConfParS(TString name) { return conf_strings[name]; }
    void InitEvent(SNDLHCEventHeader* e)
    {
        eventHeader = e;   // get mapping to eventHeader
    }

    /**      Initialization of the detector is done here    */
    virtual void Initialize();

    /**       this method is called for each step during simulation
     *       (see FairMCApplication::Stepping())
     */
    virtual Bool_t ProcessHits(FairVolume* v = 0);

    /**       Registers the produced collections in FAIRRootManager.     */
    virtual void Register();

    /** Gets the produced collections */
    virtual TClonesArray* GetCollection(Int_t iColl) const;

    /**      has to be called after each event to reset the containers      */
    virtual void Reset();

    /**      This method is an example of how to add your own point
     *       of type muonPoint to the clones array
     */
    ScifiPoint* AddHit(Int_t trackID,
                       Int_t detID,
                       TVector3 pos,
                       TVector3 mom,
                       Double_t time,
                       Double_t length,
                       Double_t eLoss,
                       Int_t pdgCode);

    /** The following methods can be implemented if you need to make
     *  any optional action in your detector during the transport.
     */

    virtual void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset) { ; }
    virtual void SetSpecialPhysicsCuts() { ; }
    virtual void EndOfEvent();
    virtual void FinishPrimary() { ; }
    virtual void FinishRun() { ; }
    virtual void BeginPrimary() { ; }
    virtual void PostTrack() { ; }
    virtual void PreTrack() { ; }
    virtual void BeginEvent() { ; }

    Scifi(const Scifi&);
    Scifi& operator=(const Scifi&);

    ClassDef(Scifi, 3)

        private
        :

        /** Track information to be stored until the track leaves the
         active volume.
         */
        Int_t fTrackID;                                                  //!  track index
    Int_t fVolumeID;                                                     //!  volume id
    TLorentzVector fPos;                                                 //!  position at entrance
    TLorentzVector fMom;                                                 //!  momentum at entrance
    Double32_t fTime;                                                    //!  time
    Double32_t fLength;                                                  //!  length
    Double32_t fELoss;                                                   //!  energy loss
    std::map<Int_t, std::map<Int_t, std::array<float, 2>>> fibresSiPM;   //! mapping of fibres to SiPM channels
    std::map<Int_t, std::map<Int_t, std::array<float, 2>>> siPMFibres;   //! inverse mapping
    std::map<Int_t, float> SiPMPos;                                      //! local SiPM channel position
    /** container for data points */
    TClonesArray* fScifiPointCollection;
    /** configuration parameters **/
    std::map<TString, Float_t> conf_floats;
    std::map<TString, Int_t> conf_ints;
    std::map<TString, TString> conf_strings;
    SNDLHCEventHeader* eventHeader;

  protected:
    Int_t InitMedium(const char* name);
};

#endif
