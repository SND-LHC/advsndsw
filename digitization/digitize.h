#ifndef ADVSND_DIGITIZE_H
#define ADVSND_DIGITIZE_H

#include "TGeoNavigator.h"
#include "AdvTargetPoint.h"
#include "AdvTargetHit.h"
#include "AdvMuFilterPoint.h"
#include "AdvMuFilterHit.h"
#include "Hit2MCPoints.h"
#include "TString.h"

#include <vector>
#include <map>
#include <algorithm>
#include <type_traits>
#include <iostream>

namespace advsnd {

    template <typename T>
    int GetStripId(const T& point, TGeoNavigator* nav) {
        auto detID = point.GetDetectorID();
        TString path;

        if constexpr (std::is_same<T, AdvTargetPoint>::value) {
            path = TString::Format(
            "/cave_1/Detector_0/volAdvTarget_0/Target_Layer_%d/SensorModule_%d/Target_SensorVolume_%d",
            point.GetLayer(), point.GetModule(), detID);
        }
        else if constexpr (std::is_same<T, AdvMuFilterPoint>::value) {
            path = TString::Format(
            "/cave_1/Detector_0/volAdvMuFilter_0/HCAL_Layer_%d/SensorModule_%d/HCAL_SensorVolume_%d",
            point.GetLayer(), point.GetModule(), detID);
        }
        else {
            std::cerr << "Invalid type\n";
        }

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

    template <typename T, typename U>
    std::vector<U> Digitize(const std::vector<T>& v_points, TGeoNavigator* nav) {
        std::map<int, std::vector<const T*>> hit_collector{};
        std::vector<U> v_hits;

        // --- Loop over all target points ---
        for (const auto& point : v_points) {
            int strip = GetStripId(point, nav);
            auto detector_id = point.GetDetectorID() - 999 + strip;
            hit_collector[detector_id].emplace_back(&point);
        }
        // --- Build hits ---
        for (const auto& [detector_id, points] : hit_collector) {
            v_hits.emplace_back(U(detector_id, points));
        }
        return v_hits;
    }

    template <typename T>
    Hit2MCPoints McLink(const std::vector<T>& v_points, TGeoNavigator* nav) {
        int point_index{0};
        std::vector<int> v_detector_id;
        Hit2MCPoints mc_links;
        std::map<int, std::map<int, double>> mc_points{};
        std::map<int, double> norm{};

        // --- Loop over all points ---
        for (const auto& point : v_points) {
            int strip = GetStripId(point, nav);
            const int detector_id = point.GetDetectorID() - 999 + strip;
            v_detector_id.emplace_back(detector_id);
            mc_points[detector_id][point_index++] = point.GetEnergyLoss();
            norm[detector_id] += point.GetEnergyLoss();
        }

        // --- Build MC links ---
        for (const auto detector_id : v_detector_id) {
            for (const auto& [point_id, energy_loss] : mc_points[detector_id]) {
                mc_links.Add(detector_id, point_id, energy_loss / norm[detector_id]);
            }
        }
        return mc_links;
    }

    template <typename T, typename U>
    class DigitizePoints {
        public:
            DigitizePoints(TGeoNavigator* nav) : m_nav(nav) {}
            std::vector<U> operator()(const std::vector<T>& v_points) const { return Digitize<T, U>(v_points, m_nav); }
        private:
            TGeoNavigator* m_nav;
    };

    template <typename T>
    class LinkPointsToDigi {
        public:
            LinkPointsToDigi(TGeoNavigator* nav) : m_nav(nav) {}
            Hit2MCPoints operator()(const std::vector<T>& v_points) const { return McLink<T>(v_points, m_nav); }
        private:
            TGeoNavigator* m_nav;
    };
}

#endif