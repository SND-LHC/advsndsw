#ifndef SHIPLHC_ADVTARGETPOINT_H_
#define SHIPLHC_ADVTARGETPOINT_H_ 1

#include "FairMCPoint.h"
#include "TObject.h"
#include "TVector3.h"

class AdvTargetPoint : public FairMCPoint
{

  public:
    /** Default constructor **/
    AdvTargetPoint();

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param pos      Ccoordinates at entrance to active volume [cm]
     *@param mom      Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     **/

    AdvTargetPoint(Int_t trackID,
                   Int_t detID,
                   TVector3 pos,
                   TVector3 mom,
                   Double_t tof,
                   Double_t length,
                   Double_t eLoss,
                   Int_t pdgCode);

    /** Destructor **/
    virtual ~AdvTargetPoint();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    Int_t PdgCode() const { return fPdgCode; }
    int constexpr GetStation() { return fDetectorID >> 17; }
    int constexpr GetPlane() { return (fDetectorID >> 16) % 2; }   // 0 is X-plane, 1 is Y-pane
    int constexpr GetRow() { return (fDetectorID >> 13) % 8; }
    int constexpr GetColumn() { return (fDetectorID >> 11) % 4; }
    int constexpr GetSensor() { return (fDetectorID >> 10) % 2; }
    int constexpr GetStrip() { return (fDetectorID) % 1024; }
    int constexpr GetModule() { return advsnd::target::columns * GetRow() + 1 + GetColumn(); }
    bool constexpr isVertical() { return GetPlane() == 0; };

  private:
    Int_t fPdgCode;

    /** Copy constructor **/
    AdvTargetPoint(const AdvTargetPoint& point);
    AdvTargetPoint operator=(const AdvTargetPoint& point);
    ClassDef(AdvTargetPoint, 1)
};

#endif   // SHIPLHC_ADVTARGETPOINT_H_
