#ifndef SHIPLHC_ADVTARGETHIT_H_
#define SHIPLHC_ADVTARGETHIT_H_

#include "SiSensor.h"
#include "SndlhcHit.h"
#include "TArrayD.h"

class AdvTargetPoint;

class TArrayD;

class AdvTargetHit : public SndlhcHit
{
  public:
    /** Default constructor **/
    AdvTargetHit();
    explicit AdvTargetHit(Int_t detID);

    // Constructor from AdvTargetPoint
    AdvTargetHit(Int_t detID, const std::vector<AdvTargetPoint*>&);

    /** Destructor **/
    ~AdvTargetHit() = default;

    /** Output to screen **/
    void Print() const;
    bool isValid() const { return flag; }
    bool isMasked(Int_t i) const { return fMasked[i]; }
    void SetMasked(Int_t i) { fMasked[i] = kTRUE; }
    int constexpr GetStation() { return fDetectorID >> 17; }
    int constexpr GetPlane() { return (fDetectorID >> 16) % 2; }   // 0 is X-plane, 1 is Y-pane
    int constexpr GetRow() { return (fDetectorID >> 13) % 8; }
    int constexpr GetColumn() { return (fDetectorID >> 11) % 4; }
    int constexpr GetSensor() { return (fDetectorID >> 10) % 2; }
    int constexpr GetStrip() { return (fDetectorID) % 1024; }
    int constexpr GetModule() { return advsnd::target::columns * GetRow() + 1 + GetColumn(); }
    bool constexpr isVertical() { return GetPlane() == 0; };

  private:
    bool flag;          ///< flag
    bool fMasked[16];   /// masked signal
    std::vector<Double_t> EFluct;
    int EFluctSize;

    ClassDef(AdvTargetHit, 1);
};

#endif   // SHIPLHC_ADVTARGETHIT_H_
