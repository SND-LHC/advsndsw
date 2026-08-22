#ifndef SNDHLLHC_SISTRIP_RAWTODIGI_H
#define SNDHLLHC_SISTRIP_RAWTODIGI_H

#include "SiStripIO.h"
#include "AdvHit.h"
#include "SiStripHardwareConstants.h"
#include "SiStripFEDBuffer.h"
#include "SiStripFEDChannel.h"
#include "SiStripDetInfo.h"
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <iostream>

template <uint8_t num_words>
constexpr uint16_t getADC_W(const uint8_t* data, uint_fast16_t offset, uint8_t bits_shift) {
    // get ADC from one or two bytes (at most 10 bits), and shift if needed
    return (data[offset ^ 7] + (num_words == 2 ? ((data[(offset + 1) ^ 7] & 0x03) << 8) : 0)) << bits_shift;
}

class SiStripRawToDigi {
    public:
        SiStripRawToDigi(const std::string& detinfo_file_name) : detector_info_(GetDetectorInfo(detinfo_file_name)) { ; }
        std::vector<AdvHit> operator()(const edm::Wrapper<FEDRawDataCollection>& sistrip_raw) const;
    private:
        std::vector<DetectorInfo> detector_info_;
};

std::vector<AdvHit> SiStripRawToDigi::operator()(const edm::Wrapper<FEDRawDataCollection>& sistrip_raw) const {
    std::vector<AdvHit> digis;
    // Create set of unique def_ids
    std::unordered_set<size_t> fed_ids;
    std::for_each(detector_info_.begin(), detector_info_.end(), [&] (const DetectorInfo& d){ fed_ids.insert(static_cast<size_t>(d.fedid)); });

    for (const auto fed_id : fed_ids) {
        const FEDRawData& raw_data = sistrip_raw.obj.data_[fed_id];
        FEDBuffer buffer(raw_data);
        if (!buffer.isValid()) {
            std::cout << "Buffer is not valid, skipping FED " << fed_id << "\n";
            continue;
        }
        if (!buffer.isZeroSuppressed()) {
            std::cout << "Only Zero Suppressed mode is supported, skipping FED " << fed_id << "\n";
            continue;
        }

        // Map fed channel to detector info
        std::unordered_map<int, const DetectorInfo*> index;
        for (const auto& d : detector_info_) {
            if (static_cast<size_t>(d.fedid) == fed_id) {
                auto [it, inserted] = index.emplace(d.fedchannel, &d);
                if (!inserted) {
                    throw std::runtime_error("Duplicate fedchannel in DetectorInfo list");
                }
            }
        }

        buffer.findChannels();
        // Loop on FED channels
        for (uint8_t i_ch{0}; i_ch < stripsensor::FEDCH_PER_FED; ++i_ch) {
            // TODO skip bad channels
            const auto channel = buffer.channel(i_ch);
            if (channel.length() == 0) {
                continue;
            }
            const uint32_t fed_key = ((fed_id & 0xFFFF) << 16) | (i_ch & 0xFFFF);
            auto it_detinfo = index.find(i_ch);
        
            constexpr uint16_t stripStart{0};
            constexpr uint8_t num_words{1};
            constexpr uint8_t headerLength{7};
            constexpr uint8_t bits_shift{0};
            const uint8_t* const data = channel.data();
            uint_fast16_t offset = channel.offset() + headerLength;
            uint_fast8_t firstStrip{0}, nInCluster{0}, inCluster{0};
            const uint_fast16_t end = channel.offset() + channel.length();
            while (offset != end) {
                if (inCluster == nInCluster) {
                    if (offset + 2 >= end) {
                        // offset should already be at end then (empty cluster)
                        break;
                    }
                    const uint_fast8_t newFirstStrip = data[(offset++) ^ 7];
                    firstStrip = newFirstStrip;
                    nInCluster = data[(offset++) ^ 7];
                    inCluster = 0;
                }

                if (it_detinfo == index.end()) throw std::runtime_error("fedchannel not found");
                auto& detector_info = *(it_detinfo->second);
                
                // The strip id in a module ranges 0 - 756, depending on the APV
                const uint16_t module_strip_id = stripsensor::SISTRIPS_PER_APV_PAIR * GetApvPair(detector_info) + (stripStart + firstStrip + inCluster);

                // For the moment set time to 0
                digis.emplace_back(AdvHit(module_strip_id, getADC_W<num_words>(data, offset, bits_shift), 0, detector_info));
                offset += num_words;
                ++inCluster;
            }
        }
    }
    // Do we need to order them? (probably not)
    return digis;
}

#endif