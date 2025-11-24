#ifndef SHIPLHC_ADVMUFILTERPOINT_H_
#define SHIPLHC_ADVMUFILTERPOINT_H_ 1

#include "FairMCPoint.h"
#include "SiSensor.h"
#include "TObject.h"
#include "TVector3.h"

class AdvMuFilterPoint : public FairMCPoint
{

  public:
    /** Default constructor **/
    AdvMuFilterPoint();

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param pos      Ccoordinates at entrance to active volume [cm]
     *@param mom      Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     **/

    AdvMuFilterPoint(Int_t trackID,
                     Int_t detID,
                     TVector3 pos,
                     TVector3 mom,
                     Double_t tof,
                     Double_t length,
                     Double_t eLoss,
                     Int_t pdgCode);

    /** Destructor **/
    virtual ~AdvMuFilterPoint();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    Int_t PdgCode() const { return fPdgCode; }
    int constexpr GetLayer() const { return fDetectorID >> 17; }
    int constexpr GetPlane() const { return (fDetectorID >> 16) % 2; }   // 0 is X-plane, 1 is Y-pane
    int constexpr GetRow() const { return (fDetectorID >> 13) % 8; }
    int constexpr GetColumn() const { return (fDetectorID >> 11) % 4; }
    int constexpr GetSensor() const { return (fDetectorID >> 10) % 2; }
    int constexpr GetStrip() const { return (fDetectorID) % 1024; }
    int constexpr GetModule() const { return advsnd::hcal::columns * GetRow() + 1 + GetColumn(); }
    bool constexpr isVertical() const { return GetPlane() == 1; };

  private:
    Int_t fPdgCode;

    ClassDef(AdvMuFilterPoint, 1)
};

#endif   // SHIPLHC_ADVMUFILTERPOINT_H_
