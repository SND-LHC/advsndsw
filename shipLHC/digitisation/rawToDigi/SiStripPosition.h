#ifndef SNDHLLHC_SISTRIP_POSITION_H
#define SNDHLLHC_SISTRIP_POSITION_H

#include <cmath>
#include <thread>
#include "Math/Point3D.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"

#include "SiStripHardwareConstants.h"

ROOT::Math::XYZPoint GetSiStripPosition(uint32_t detector_id) {

    const int strip = detector_id & 0x3FF;                // actual strip ID
    const int geofile_detID = detector_id - strip + 999;   // the det id number needed to read the geometry
 
    const int layer = (geofile_detID >> 13);
    const int row = (geofile_detID >> 11) & 0x3;
    const int column = (geofile_detID >> 10) & 0x1;
    const int module_id = geofile_detID;

    TString path = TString::Format("/cave_1/"
                                   "Detector_0/"
                                   "volAdvTarget_0/"
                                   "Target_Layer_%d/"
                                   "Row_%d_Column_%d_0/"
                                   "Target_DoubleSensorVolume_%d",
                                   layer,
                                   row,
                                   column,
                                   module_id);
                                   
    thread_local TGeoNavigator* nav = gGeoManager->AddNavigator();

    if (!nav) {
        std::cerr << "NULL navigator in thread " << std::this_thread::get_id() << std::endl;
        return {std::nan(""), std::nan(""), std::nan("")};
    }
    if (!nav->cd(path)) {
        std::cerr<<"Wrong path " << path << "\n";
        return {std::nan(""), std::nan(""), std::nan("")};
    }

    // Knowing the strip, get the postion along the module
    const double strip_offset = static_cast<double>(strip) - MAX_SISTRIPS_PER_MODULE / 2.0;
    double local_mid_pos[3] = {0.0, strip_offset * SENSOR_LENGTH_CM / MAX_SISTRIPS_PER_MODULE, 0.0};
    double global_mid_pos[3];
    nav->LocalToMaster(local_mid_pos, global_mid_pos);

    // Return SiStrip midpoint
    return {global_mid_pos[0], global_mid_pos[1], global_mid_pos[2]};
}

#endif