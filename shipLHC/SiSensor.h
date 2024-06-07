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
        const int rows = 8;
        const int columns = 4;
        const int planes = 2;
        const double module_row_gap = 0.5 * ShipUnit::mm;
        const double module_column_gap = 13.9 * ShipUnit::mm;
    }

    namespace muon {
        const int rows = 8;
        const int columns = 4;
        const int planes = 1;
        const double module_row_gap = 0.5 * ShipUnit::mm;
        const double module_column_gap = 13.9 * ShipUnit::mm;
    }

}

#endif // SISENSOR_H_
