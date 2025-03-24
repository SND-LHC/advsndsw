#ifndef SHIPLHC_ADVTARGETPOINT_H_
#define SHIPLHC_ADVTARGETPOINT_H_ 1

#include "FairMCPoint.h"
#include "SiSensor.h"
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
     *@param pos      Coordinates halfway between entrance and exit of active volume [cm]
     *@param mom      Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     *@param pdgCode  PDG code of MC particle
     *@param exitpoint      Coordinates at exit from active volume [cm]
     **/

    AdvTargetPoint(Int_t trackID,
                   Int_t detID,
                   TVector3 pos,
                   TVector3 mom,
                   Double_t tof,
                   Double_t length,
                   Double_t eLoss,
                   Int_t pdg_code,
                   TVector3 exit_point);

    /** Destructor **/
    virtual ~AdvTargetPoint();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    Int_t PdgCode() const { return fPdgCode; }
    int constexpr GetLayer() { return fDetectorID >> 17; }
    int constexpr GetPlane() { return (fDetectorID >> 16) % 2; }   // 0 is X-plane, 1 is Y-pane
    int constexpr GetRow() { return (fDetectorID >> 13) % 8; }
    int constexpr GetColumn() { return (fDetectorID >> 11) % 4; }
    int constexpr GetSensor() { return (fDetectorID >> 10) % 2; }
    int constexpr GetStrip() { return (fDetectorID) % 1024; }
    int constexpr GetModule() { return advsnd::target::columns * GetRow() + 1 + GetColumn(); }
    bool constexpr IsVertical() { return GetPlane() == 1; };
    TVector3 GetEntryPoint() const { return TVector3(2 * fX - fExitX, 2 * fY - fExitY, 2 * fZ - fExitZ); }
    TVector3 GetExitPoint() const { return TVector3(fExitX, fExitY, fExitZ); }

  private:
    Int_t fPdgCode;
    Double_t fExitX;
    Double_t fExitY;
    Double_t fExitZ;

    /** Copy constructor **/
    AdvTargetPoint(const AdvTargetPoint& point);
    AdvTargetPoint operator=(const AdvTargetPoint& point);
    ClassDef(AdvTargetPoint, 1)
};

#endif   // SHIPLHC_ADVTARGETPOINT_H_
