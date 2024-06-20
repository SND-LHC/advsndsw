#ifndef ADVTARGETHIT_H
#define ADVTARGETHIT_H 1

#include "SndlhcHit.h"

#include <map>

class AdvTargetPoint;

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
    int constexpr GetStation() { return floor(fDetectorID / 1e7); }
    int constexpr GetPlane() { return int(fDetectorID / 1e6) % 2; }   // 0 is X-plane, 1 is Y-pane
    int constexpr GetRow() { return int(fDetectorID / 1e5) % 10; }
    int constexpr GetColumn() { return int(fDetectorID / 1e5) % 10; }
    int constexpr GetSensor() { return int(fDetectorID / 1e4) % 10; }
    int constexpr GetStrip() { return int(fDetectorID) % 1000; }
    int constexpr GetModule() { return 2 * GetRow() + 1 + GetColumn(); }
    bool constexpr isVertical() { return GetPlane() == 0; };

  private:
    bool flag;          ///< flag
    bool fMasked[16];   /// masked signal

    ClassDef(AdvTargetHit, 1);
};

#endif
