#ifndef SHIPLHC_ADVPOINT_H_
#define SHIPLHC_ADVPOINT_H_ 1

#include "FairMCPoint.h"
#include "SiSensor.h"
#include "TObject.h"
#include "TVector3.h"

class AdvPoint : public FairMCPoint
{

  public:
    /** Default constructor **/
    AdvPoint();

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

    AdvPoint(Int_t trackID,
                   Int_t detID,
                   TVector3 pos,
                   TVector3 mom,
                   Double_t tof,
                   Double_t length,
                   Double_t eLoss,
                   Int_t pdg_code,
                   TVector3 exit_point);

    /** Destructor **/
    virtual ~AdvPoint();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    Int_t PdgCode() const { return fPdgCode; }
    int constexpr GetLayer() const { return fDetectorID >> 13; }
    int constexpr GetRow() const { return (fDetectorID >> 11) & 0x3; }
    int constexpr GetColumn() const { return (fDetectorID >> 10) & 0x1; }
    int constexpr GetStrip() const { return (fDetectorID) & 0x3FF; }
    int constexpr GetModule(int system, int setup = 0) const
    {
        if (system == 1)
        {
          return (setup == 0) ? (advsnd::target::columns * GetRow() + 1 + GetColumn())
                              : (advsnd::tb_target::columns * GetRow() + 1 + GetColumn());
        }
        else
        {
          return advsnd::hcal::columns * GetRow() + 1 + GetColumn();
        }
    }
    bool constexpr IsVertical() const { return GetLayer() % 2 == 0; }; // 1 is X-plane, 0 is Y-pane for TB26
    TVector3 GetEntryPoint() const { return TVector3(2 * fX - fExitX, 2 * fY - fExitY, 2 * fZ - fExitZ); }
    TVector3 GetExitPoint() const { return TVector3(fExitX, fExitY, fExitZ); }

  private:
    Int_t fPdgCode;
    Double_t fExitX;
    Double_t fExitY;
    Double_t fExitZ;

    ClassDef(AdvPoint, 1)
};

#endif   // SHIPLHC_ADVPOINT_H_
