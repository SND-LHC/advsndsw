#ifndef SHIPLHC_ADVHIT_H_
#define SHIPLHC_ADVHIT_H_

#include "SiSensor.h"
#include "TObject.h"
#include <vector>
#include <cstdint>

class AdvHit : public TObject
{
  public:
    /** Default constructor **/
    AdvHit() : detector_id_(0), daq_id_(0), time_(0), signal_(0.0), is_valid_(true) {}

    explicit AdvHit(uint32_t detID) : detector_id_(detID), daq_id_(0), time_(0), signal_(0.0), is_valid_(true) {}

    /** Destructor **/
    ~AdvHit() = default;

    /** Output to screen **/
    void Print() const;
    /** Setters **/
    void SetSignal(uint16_t set_signal) { signal_ = set_signal; }
    void SetTime(float set_time) { time_ = set_time; }
    void SetValid(bool is_valid) { is_valid_ = is_valid; }
    /** Getters **/
    uint16_t GetSignal() const { return signal_; }
    float GetTime() const { return time_; }
    int GetDetectorID() const { return detector_id_; }
    bool IsValid() const { return is_valid_; }
    int constexpr GetLayer() const { return detector_id_ >> 13; }
    int constexpr GetRow() const { return (detector_id_ >> 11) & 0x3; }
    int constexpr GetColumn() const { return (detector_id_ >> 10) & 0x1; }
    int constexpr GetStrip() const { return (detector_id_) & 0x3FF; }
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

  private:
    uint32_t detector_id_;
    uint32_t daq_id_;
    float time_;
    uint16_t signal_;
    bool is_valid_;

    ClassDef(AdvHit, 2);
};

#endif   // SHIPLHC_ADVHIT_H_
