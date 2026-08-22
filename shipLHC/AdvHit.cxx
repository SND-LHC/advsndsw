#include "AdvHit.h"
#include "SiStripDetInfo.h"

#include <iostream>

// -----   Public method Print   -------------------------------------------
void AdvHit::Print() const
{
    std::cout << " AdvHit: in detector " << detector_id_ << "\tsignal: " << signal_ << " \ttime: " << time_ << "\n";
}

// Constructor from raw data, currntly no time info is retrieved
AdvHit::AdvHit(uint16_t strip, uint16_t adc, float time, const DetectorInfo& detinfo) : daq_id_(detinfo.ccuaddress), time_(time), signal_(adc), is_valid_(true) {
  detector_id_ =
    ((detinfo.layer & 0x7FFFF) << 13) |
    ((detinfo.row & 0x3) << 11) |
    ((detinfo.column & 0x1) << 10) |
    (strip & 0x3FF);
}