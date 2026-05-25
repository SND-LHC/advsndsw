#ifndef SHIPLHC_ADVHIT_H_
#define SHIPLHC_ADVHIT_H_

#include "TObject.h"

#include "SiSensor.h"
#include "digitisation/AdvSignal.h"
#include <map>

#ifndef __CINT__
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#endif //__CINT__

class AdvPoint;

class AdvHit : public TObject
{
  public:
    /** Default constructor **/
    AdvHit();
    explicit AdvHit(Int_t detID);

    // Constructor from AdvPoint
    AdvHit(Int_t detID, const std::vector<AdvPoint*>&);

    /** Destructor **/
    ~AdvHit() = default;

    /** Output to screen **/
    void Print() const;
    //bool isValid() const { return flag; }
    /** Setters **/
    void SetSignal(float set_signal) { signal=set_signal; }
    void SetTime(float set_time) { time = set_time; }
    /** Getters **/
    float GetSignal() { return signal; }
    float GetTime() { return time; }
    int GetDetectorID() { return fDetectorID; }
    //bool isMasked(Int_t i) const { return fMasked[i]; }
    //void SetMasked(Int_t i) { fMasked[i] = kTRUE; }
    std::map<std::string, std::vector<Int_t>> GetHit() { return fDigitisedHit; }
    int constexpr GetLayer() { return fDetectorID >> 13; }
    int constexpr GetRow() { return (fDetectorID >> 11) & 0x3; }
    int constexpr GetColumn() { return (fDetectorID >> 10) & 0x1; }
    int constexpr GetStrip() { return (fDetectorID) & 0x3FF; }
    int constexpr GetModule(int system, int setup = 0)
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
    bool constexpr IsVertical() { return GetLayer() % 2 == 0; }; // 1 is X-plane, 0 is Y-pane for TB26

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& boost::serialization::base_object<TObject>(*this);
        ar& fDetectorID;
        ar& fDaqID;
        ar& signal;
        ar& time;
    }

  protected:
#ifndef __CINT__ // for BOOST serialization
    friend class boost::serialization::access;
#endif // for BOOST serialization

  private:
    int fDetectorID;
    int fDaqID;
    float signal;
    float time;
    bool flag;          ///< flag
    std::map<std::string, std::vector<Int_t>> fDigitisedHit;  //! don't store this item
    ClassDef(AdvHit, 1);
};

#endif   // SHIPLHC_ADVHIT_H_
