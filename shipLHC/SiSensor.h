#ifndef SISENSOR_H_
#define SISENSOR_H_

/* Header-only to implement common code for CMS TOB */
/* #include "TGeoBBox.h" */
/* #include "TGeoVolume.h" */
#include "ShipUnit.h"

namespace advsnd {

    const double sensor_width = 93.7 * ShipUnit::mm;
    const double sensor_length = 91.5 * ShipUnit::mm;
    const double module_length = 23.95 * ShipUnit::cm;
    const double module_width = 12.0 * ShipUnit::cm;
    const double sensor_gap = 3.1 * ShipUnit::mm;
    const int sensors = 2;
    const int strips = 768;

    namespace target {
        const int rows = 4;
        const int columns = 2;
        const int planes = 2;
        const double module_row_gap = 0.5 * ShipUnit::mm;
        const double module_column_gap = 13.9 * ShipUnit::mm;
    }

    namespace hcal {
        const int rows = 6;
        const int columns = 3;
        const int planes = 2;
        const int stations = 5;
        const double module_row_gap = -26.3 * ShipUnit::mm + sensor_gap;
        const double module_column_gap = 2 * (sensor_length + sensor_gap) - module_length;
        const double plane_width = (2 * columns) * sensor_length + (2 * columns - 1) * sensor_gap + 2 * 6.45 * ShipUnit::mm;
    }

    namespace muon {
        const int rows = 6;
        const int columns = 3;
        const int planes = 1;
        const double module_row_gap = -26.3 * ShipUnit::mm + sensor_gap;
        const double module_column_gap = 2 * (sensor_length + sensor_gap) - module_length;
        const double plane_width = (2 * columns) * sensor_length + (2 * columns - 1) * sensor_gap + 2 * 6.45 * ShipUnit::mm;
    }

}

#endif // SISENSOR_H_
