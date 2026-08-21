#ifndef SNDHLLHC_SISTRIP_FEDFULLDEBUGHEADER_H
#define SNDHLLHC_SISTRIP_FEDFULLDEBUGHEADER_H

#include "SiStripHardwareConstants.h"

#include <cstdint>
#include <cstring>
#include <memory>

class FEDFullDebugHeader {
  public:
    FEDFullDebugHeader(const uint8_t* headerBuffer) { memcpy(header_, headerBuffer, FULL_DEBUG_HEADER_SIZE_IN_BYTES);}
    inline uint16_t feUnitLength(const uint8_t internalFEUnitNum) const { return ((feWord(internalFEUnitNum)[15] << 8) | (feWord(internalFEUnitNum)[14]));}
    inline bool fePresent(const uint8_t internalFEUnitNum) const { return (feUnitLength(internalFEUnitNum) != 0);}
    inline const uint8_t* feWord(const uint8_t internalFEUnitNum) const { return header_ + internalFEUnitNum * 2 * 8;}
    inline uint32_t daqRegister() const { return get32BitWordFrom(feWord(7) + 10); }
    inline static uint32_t get32BitWordFrom(const uint8_t* startOfWord) { return (startOfWord[0] | (startOfWord[1] << 8) | (startOfWord[2] << 16) | (startOfWord[3] << 24));}
    static constexpr size_t FULL_DEBUG_HEADER_SIZE_IN_64BIT_WORDS = FEUNITS_PER_FED * 2;
    static constexpr size_t FULL_DEBUG_HEADER_SIZE_IN_BYTES = FULL_DEBUG_HEADER_SIZE_IN_64BIT_WORDS * 8;
  private:
    uint8_t header_[FULL_DEBUG_HEADER_SIZE_IN_BYTES];
};
#endif