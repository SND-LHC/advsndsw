#ifndef SNDHLLHC_SISTRIP_DIGI_CLUSTERING_H
#define SNDHLLHC_SISTRIP_DIGI_CLUSTERING_H

#include "SiStripIO.h"
#include "AdvHit.h"
#include "SiStripHardwareConstants.h"
#include <vector>
#include <map>
#include <cstdint>

class SiStripDigiClustering {
    public:
        SiStripDigiClustering() { ; }
        SiStripClusteringProducts operator()(const std::vector<AdvHit>& input_digis) const;
};

SiStripClusteringProducts SiStripDigiClustering::operator()(const std::vector<AdvHit>& input_digis) const
{
    std::map<Module, std::vector<AdvHit>> grouped;
    std::vector<SiStripCluster> clusters;
    std::vector<AdvHit> clustered_output;
    clustered_output.reserve(input_digis.size());

    // Group by module
    for (const auto& d : input_digis)
    {
        grouped[{d.GetLayer(), d.GetRow(), d.GetColumn()}].push_back(d);
    }

    constexpr uint16_t seed_thr  = 3 * SISTRIP_NOISE_ADC;
    constexpr uint16_t neigh_thr = 2 * SISTRIP_NOISE_ADC;
    constexpr uint16_t cluster_cut_factor = SISTRIP_NOISE_ADC * SISTRIP_NOISE_ADC;

    // Cluster each module independently
    for (auto& [mod, digis] : grouped)
    {
        std::sort(digis.begin(), digis.end(), [](const AdvHit& a, const AdvHit& b) { return a.GetStrip() < b.GetStrip();});
        std::vector<uint8_t> used(digis.size(), 0);

        for (size_t i = 0; i < digis.size(); ++i)
        {
            if (used[i]) continue;

            const AdvHit& seed = digis[i];

            if (seed.GetSignal() <= seed_thr)
                continue;

            std::vector<size_t> cluster;
            uint32_t sum_signal{0};

            // start cluster
            cluster.push_back(i);
            used[i] = 1;
            sum_signal += seed.GetSignal();

            int last_strip = seed.GetStrip();

            // forward
            for (size_t j = i + 1; j < digis.size(); ++j)
            {
                if (used[j]) continue;

                const AdvHit& d = digis[j];

                if (d.GetStrip() != last_strip + 1)
                    break;

                if (d.GetSignal() <= neigh_thr)
                    break;

                cluster.push_back(j);
                used[j] = 1;
                sum_signal += d.GetSignal();
                last_strip = d.GetStrip();
            }

            // backward
            last_strip = seed.GetStrip();

            for (int j = (int)i - 1; j >= 0; --j)
            {
                if (used[j]) continue;

                const AdvHit& d = digis[j];

                if (d.GetStrip() != last_strip - 1)
                    break;

                if (d.GetSignal() <= neigh_thr)
                    break;

                cluster.push_back(j);
                used[j] = 1;
                sum_signal += d.GetSignal();
                last_strip = d.GetStrip();
            }

            uint32_t threshold = static_cast<uint32_t>(cluster.size()) * cluster_cut_factor;

            if (sum_signal > threshold)
            {   
                size_t digi_with_most_signal{cluster.front()};
                for (auto idx : cluster) {
                    clustered_output.push_back(digis[idx]);
                    if (digis[idx].GetSignal() > digis[digi_with_most_signal].GetSignal()) {
                        digi_with_most_signal = idx;
                    }
                }
                clusters.emplace_back(digis[digi_with_most_signal].GetDetectorID(), sum_signal, cluster.size(), mod.layer, mod.row, mod.col, digis[digi_with_most_signal].IsVertical());
            }
            else
            {
                // reject cluster → free digis
                for (auto idx : cluster)
                    used[idx] = 0;
            }
        }
    }
    return {std::move(clustered_output), std::move(clusters)};
}

#endif