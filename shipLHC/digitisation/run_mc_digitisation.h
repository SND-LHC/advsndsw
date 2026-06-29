#ifndef ADVSND_RUN_MC_DIGITISATION_H
#define ADVSND_RUN_MC_DIGITISATION_H

#include "TGeoNavigator.h"
#include "AdvPoint.h"
#include "AdvHit.h"
#include "Hit2MCPoints.h"
#include "TString.h"
#include "AdvDigitisation.h"
#include "SiDigiParameters.h"

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

namespace advsnd {

    // Adds strip saturation
    int Saturate(int signal)
    {
        if (stripsensor::frontend::ZSModeOption)
        {
            if (signal > 1022) signal = 255;
            else if (signal > 253) signal = 254;
            else if (signal < 0) signal = 0;
        } else {
            if (signal > 1023) signal = 1023;
            else if (signal < 0) signal = 0;
        }
        return signal;
    }

    int GetStripId(const AdvPoint& point, TGeoNavigator* nav) {
        auto path = TString::Format("/cave_1/"
                                    "Detector_0/"
                                    "volAdvTarget_0/"
                                    "Target_Layer_%d/"
                                    "Row_%d_Column_%d_0/"
                                    "Target_DoubleSensorVolume_%d",
                                    point.GetLayer(),
                                    point.GetRow(),
                                    point.GetColumn(),
                                    point.GetDetectorID());

        if (nav->CheckPath(path)) {
            nav->cd(path);
        } else {
            std::cerr << "Invalid geometry path: " << path << "\n";
        }

        nav->GetCurrentNode();
        double global_pos[3] = {point.GetX(), point.GetY(), point.GetZ()};
        double local_pos[3];
        nav->MasterToLocal(global_pos, local_pos);

        int strip = floor((local_pos[1] / (advsnd::sensor_length / advsnd::strips)) + (advsnd::strips / 2));
        strip = std::clamp(strip, 0, advsnd::strips - 1);
        return strip;
    }

    std::vector<AdvHit> Digitize(const std::vector<AdvPoint>& v_points, TGeoNavigator* nav) {
        std::map<int, std::vector<const AdvPoint*>> hit_collector{};
        std::vector<AdvHit> v_hits;
        std::map<int, int> module_map{};
        std::map<int, std::vector<std::map<std::string, std::vector<Int_t>>>> module_collector{};

        AdvDigitisation advdigi{};

        for (const auto& point : v_points) {
            int strip = GetStripId(point, nav);
            auto detector_id = point.GetDetectorID() - 999 + strip;

            // make a list of detector_id in the same unique module
            module_map[detector_id] = point.GetDetectorID();
            // Collect points by virtual strip
            hit_collector[detector_id].emplace_back(&point);
        }
        for (const auto& [detector_id, points] : hit_collector) {
            // Make one hit per virtual strip (detector ID module + strip)
            module_collector[module_map[detector_id]].emplace_back(advdigi.digirunoutput(detector_id, points));
        }

        // Loop over hits in an unique module and sum the charges per unique strip
        // to form its total integrated charge
        for (const auto& [detID, digihits] : module_collector) {
            std::vector<float> sum_adc(advsnd::strips, 0);
            for (const auto& fDigitisedHit : digihits) {
                for(int a =0; a<fDigitisedHit.at("Strips").size(); a++) {
                    sum_adc[fDigitisedHit.at("Strips")[a]] += fDigitisedHit.at("ADC")[a];
                }
            }// end loop over hits in the same module.
            // At this stage one has the total changer per strip in a module.
            // Now one writes the ADC to the respective digi hit, respecting saturation!
            std::vector<int> existing_hit{};
            for (const auto& fDigitisedHit : digihits) {
                auto it = std::find_if(std::begin(module_map), std::end(module_map),
                    [&detID]( const auto &p )
                    {
                        return p.second == detID; 
                    } );
                    
                if (it == std::end(module_map)) 
                {
                    std::cout << "detector_id not found, skipping hit.\n";
                    continue;
                }
                int detector_id = it->first;
                int strip = detector_id & 0x3FF;
                existing_hit.push_back(strip);
            }
            // Add new hits for all strips that have "non-zero" charge, but were not intersected by a particle
            // meaning one now adds noise hits.
            // Beware "non-zero charge" in Virgin Raw is real non-zero, while for Zero Suppression it means charge
            // below a pre-set threshold.
            for (auto i=0; i < advsnd::strips; i++)
            {
                if (sum_adc[i] == 0) continue;
                if (count(existing_hit.begin(), existing_hit.end(), i)) continue;
                else {
                    auto detector_id = detID - 999 + i;
                    v_hits.emplace_back(AdvHit(detector_id));
                    v_hits.front().SetSignal(Saturate(sum_adc[i]));
                }
            }
        }
        return v_hits;
    }

    Hit2MCPoints McLink(const std::vector<AdvPoint>& v_points, TGeoNavigator* nav) {
        int point_index{0};
        Hit2MCPoints mc_links;
        std::vector<int> v_detector_id;
        std::map<int, std::map<int, double>> mc_points{};
        std::map<int, double> norm{};

        for (const auto& point : v_points) {
            int strip = GetStripId(point, nav);
            const int detector_id = point.GetDetectorID() - 999 + strip;

            v_detector_id.emplace_back(detector_id);
            mc_points[detector_id][point_index++] = point.GetEnergyLoss();
            norm[detector_id] += point.GetEnergyLoss();
        }

        for (const auto detector_id : v_detector_id) {
            auto point_map = mc_points[detector_id];
            for (const auto& [point_id, energy_loss] : point_map) {
                mc_links.Add(detector_id, point_id, energy_loss / norm[detector_id]);
            }
        }
        return mc_links;
    }


    class DigitizePoints {
        public:
            DigitizePoints(TGeoNavigator* nav) : m_nav(nav) {}
            std::vector<AdvHit> operator()(const std::vector<AdvPoint>& v_points) const { return Digitize(v_points, m_nav); }
        private:
            TGeoNavigator* m_nav;
    };

    class LinkPointsToDigi {
        public:
            LinkPointsToDigi(TGeoNavigator* nav) : m_nav(nav) {}
            Hit2MCPoints operator()(const std::vector<AdvPoint>& v_points) const { return McLink(v_points, m_nav); }
        private:
            TGeoNavigator* m_nav;
    };
}

#endif