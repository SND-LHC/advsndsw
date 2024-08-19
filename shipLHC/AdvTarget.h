//
//  AdvTarget.h
//
//  O. Lantwin and D. Centanni
//  Dec 2023
//
#ifndef SHIPLHC_ADVTARGET_H_
#define SHIPLHC_ADVTARGET_H_

#include "FairDetector.h"
#include "Rtypes.h"   // for ShipMuonShield::Class, Bool_t, etc
#include "TLorentzVector.h"

class AdvTargetPoint;
class FairVolume;
class TVector3;
class TClonesArray;

class AdvTarget : public FairDetector
{
  public:
    AdvTarget(const char *name, Bool_t Active, const char *Title = "AdvTarget");
    AdvTarget();
    virtual ~AdvTarget();

    /** Create the detector geometry **/
    void ConstructGeometry();

    void SetConfPar(TString name, Float_t value) { conf_floats[name] = value; }
    void SetConfPar(TString name, Int_t value) { conf_ints[name] = value; }
    void SetConfPar(TString name, TString value) { conf_strings[name] = value; }
    /** End-point positions of strips, centre line, associated with readout channel **/
    void GetPosition(Int_t detID, TVector3 &A, TVector3 &B);
    Float_t GetConfParF(TString name) { return conf_floats[name]; }
    Int_t GetConfParI(TString name) { return conf_ints[name]; }
    TString GetConfParS(TString name) { return conf_strings[name]; }

    /** Initialization of the detector is done here    */
    virtual void Initialize();

    /**       this method is called for each step during simulation
     *       (see FairMCApplication::Stepping())
     */
    virtual Bool_t ProcessHits(FairVolume *v = 0);

    /**       Registers the produced collections in FAIRRootManager.    */
    virtual void Register();

    /** Gets the produced collections */
    virtual TClonesArray *GetCollection(Int_t iColl) const;

    /**      has to be called after each event to reset the containers      */
    virtual void Reset();

    /**      This method adds an AdvTargetPoint to the clones array
     */
    AdvTargetPoint *AddHit(Int_t trackID,
                           Int_t detID,
                           TVector3 pos,
                           TVector3 mom,
                           Double_t time,
                           Double_t length,
                           Double_t eLoss,
                           Int_t pdgCode,
                           TVector3 exit_point);

    /** The following methods can be implemented if you need to make
     *  any optional action in your detector during the transport.
     */

    virtual void CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset) { ; }
    virtual void SetSpecialPhysicsCuts() { ; }
    virtual void EndOfEvent();
    virtual void FinishPrimary() { ; }
    virtual void FinishRun() { ; }
    virtual void BeginPrimary() { ; }
    virtual void PostTrack() { ; }
    virtual void PreTrack() { ; }
    virtual void BeginEvent() { ; }

    AdvTarget(const AdvTarget &);
    AdvTarget &operator=(const AdvTarget &);

  private:
    /** Track information to be stored until the track leaves the
        active volume.
    */
    Int_t fTrackID;               //!  track index
    Int_t fVolumeID;              //!  volume id
    TLorentzVector fEntryPoint;   //!  position at entrance
    TLorentzVector fMom;          //!  momentum at entrance
    Double32_t fTime;             //!  time
    Double32_t fLength;           //!  length
    Double32_t fELoss;            //!  energy loss
    /** container for data points */
    TClonesArray *fAdvTargetPointCollection;
    /** configuration parameters **/
    std::map<TString, Float_t> conf_floats;
    std::map<TString, Int_t> conf_ints;
    std::map<TString, TString> conf_strings;

    ClassDef(AdvTarget, 1)

        protected
        :

        Int_t InitMedium(const char *name);
};

#endif   // SHIPLHC_ADVTARGET_H_
