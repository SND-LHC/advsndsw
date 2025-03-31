#ifndef SHIPLHC_SISENSOR_H_
#define SHIPLHC_SISENSOR_H_

/* Header-only to implement common code for CMS TOB */
#include "ShipUnit.h"

namespace advsnd {

const double sensor_length = 93.7 * ShipUnit::mm;
const double sensor_width = 91.5 * ShipUnit::mm;
const double sensor_thickness = 0.5 * ShipUnit::mm;
const double module_width = 23.95 * ShipUnit::cm;
const double module_length = 12.0 * ShipUnit::cm;
const double support_thickness = 1.8 * ShipUnit::mm;
const double sensor_gap = 3.1 * ShipUnit::mm;
const int sensors = 2;
const int strips = 768;

namespace target {
const int rows = 4;
const int columns = 2;
const int planes = 2;
const double module_dead_space_side_small = 6.45 * ShipUnit::mm;
const double module_dead_space_side_large = 46.95 * ShipUnit::mm;
const double module_dead_space_top = 13.15 * ShipUnit::mm;
const double module_dead_space_bottom = module_dead_space_top;
// cenral gap between columns/rows in a detective layer
const double layer_central_gap = 13.4 * ShipUnit::mm;
// column gap btw 2 modules, excluding dead spaces
const double modules_column_gap = layer_central_gap - 2 * module_dead_space_side_small;
// 2 rows of modules are staggred to form a module 'tandem'
const double modules_rows_overlap = 27.6 * ShipUnit::mm;
// two tandems of modules are staggered to form the central gap
const double tandem_modules_rows_overlap = module_dead_space_top + module_dead_space_bottom - layer_central_gap;
const double offset = module_dead_space_side_large - module_dead_space_top;
}   // namespace target

namespace hcal {
const int rows = 4;
const int columns = 2;
const int planes = 2;
const int n_XY_layers = 28;
const int n_X_layers = 6;
const double module_dead_space_side_small = 6.45 * ShipUnit::mm;
const double module_dead_space_side_large = 46.95 * ShipUnit::mm;
const double module_dead_space_top = 13.15 * ShipUnit::mm;
const double module_dead_space_bottom = module_dead_space_top;
// cenral gap between columns in a detective layer
const double layer_central_gap = 13.4 * ShipUnit::mm;
// column gap btw 2 modules, excluding dead spaces
const double modules_column_gap = layer_central_gap - 2 * module_dead_space_side_small;
// rows of modules are staggred
const double modules_rows_overlap = module_dead_space_top + module_dead_space_bottom + modules_column_gap;   // 26.8mm
}   // namespace hcal

}   // namespace advsnd

#endif   // SHIPLHC_SISENSOR_H_
