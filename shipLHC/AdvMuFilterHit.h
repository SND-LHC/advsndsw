#ifndef ADVMUFILTERHIT_H
#define ADVMUFILTERHIT_H 1

#include "SndlhcHit.h"

#include <map>

class AdvMuFilterPoint;

class AdvMuFilterHit : public SndlhcHit
{
  public:
    /** Default constructor **/
    AdvMuFilterHit();
    explicit AdvMuFilterHit(Int_t detID);

    // Constructor from MuFilterPoint
    AdvMuFilterHit(Int_t detID, const std::vector<AdvMuFilterPoint*>&);

    /** Destructor **/
    ~AdvMuFilterHit() = default;

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
    int constexpr GetModule() { return 3 * GetRow() + 1 + GetColumn(); }
    bool constexpr isVertical() { return GetPlane() == 0; };

  private:
    bool flag;          ///< flag
    bool fMasked[16];   /// masked signal

    ClassDef(AdvMuFilterHit, 1);
};

#endif
