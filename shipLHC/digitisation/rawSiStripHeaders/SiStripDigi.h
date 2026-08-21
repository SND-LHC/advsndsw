#ifndef SNDHLLHC_SISTRIPDIGI_H
#define SNDHLLHC_SISTRIPDIGI_H

#include <cstdint>
#include "SiStripDetInfo.h"

class SiStripDigi {
public:
  SiStripDigi(uint16_t strip, uint16_t adc, uint32_t fed_key, uint32_t time, const DetectorInfo& detinfo);
  SiStripDigi() : detector_id_(0), fed_key_(0), time_(0), adc_(0), is_valid_(false) {}

  inline uint32_t GetFedKey() const { return fed_key_; }
  inline uint32_t GetDetectorId() const { return detector_id_; }
  inline int GetLayer() const { return detector_id_ >> 13; }
  inline int GetRow() const { return (detector_id_ >> 11) & 0x3; }
  inline int GetColumn() const { return (detector_id_ >> 10) & 0x1; }
  inline int GetStrip() const { return detector_id_ & 0x3FF; }
  inline uint32_t GetTime() const { return time_; }
  inline uint16_t GetSignal() const { return adc_; }
  inline bool IsVertical() const { return GetLayer() % 2 == 0; };
  inline bool IsValid() const { return is_valid_; }

  ClassDef(SiStripDigi, 1);

private:
  uint32_t detector_id_;
  uint32_t daq_id_;
  uint32_t fed_key_;
  uint32_t time_;
  uint16_t adc_;
  bool is_valid_;
};

inline SiStripDigi::SiStripDigi(uint16_t strip, uint16_t adc, uint32_t fed_key, uint32_t time, const DetectorInfo& detinfo) : daq_id_(detinfo.ccuaddress), fed_key_(fed_key), time_(time), adc_(adc), is_valid_(true) {
  detector_id_ =
    ((detinfo.layer & 0x7FFFF) << 13) |
    ((detinfo.row & 0x3) << 11) |
    ((detinfo.column & 0x1) << 10) |
    (strip & 0x3FF);
}
#endif