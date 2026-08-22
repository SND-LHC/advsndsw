#ifndef SNDHLLHC_SISTRIP_FEDBUFFER_H
#define SNDHLLHC_SISTRIP_FEDBUFFER_H

#include "SiStripIO.h"
#include "SiStripFEDChannel.h"
#include "SiStripHardwareConstants.h"
#include "FEDFullDebugHeader.h"

#include <vector>
#include <cstdint>
#include <memory>

// Class representing standard (non-spy channel) FED buffers
class FEDBuffer {
  public:
    explicit FEDBuffer(const FEDRawData& fedBuffer);
    void findChannels();
    const FEDChannel& channel(const uint8_t internalFEDChannelNum) const { return channels_[internalFEDChannelNum];}
    bool isValid() const;
    bool isZeroSuppressed() const;

  private:
    bool fePresent(uint8_t internalFEUnitNum) const;
    const uint8_t* getPointerToDataAfterTrackerSpecialHeader() const;
    const uint8_t* getPointerToByteAfterEndOfPayload() const;

    std::vector<FEDChannel> channels_;
    std::unique_ptr<FEDFullDebugHeader> feHeader_;
    const uint8_t* originalBuffer_;
    const uint8_t* orderedBuffer_;
    const uint8_t* payloadPointer_;
    size_t bufferSize_;
    uint16_t payloadLength_;
    uint8_t validChannels_;
    bool fePresent_[stripsensor::FEUNITS_PER_FED];
};

#endif