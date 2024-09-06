#ifndef SHIPLHC_SIDIGIPARAMETERS_H_
#define SHIPLHC_SIDIGIPARAMETERS_H_

/* Header-only to implement common code for CMS TOB */
#include "ShipUnit.h"

namespace stripsensor {

const bool peakmode = 1;
const std::string APVpeakpulse = "advsndsw/shipLHC/data/APVShapePeak_default.txt";
const std::string APVdecopulse = "advsndsw/shipLHC/data/APVShapeDeco_default.txt";

const Double_t value_per_mip = 1; // in mA
const Double_t amplificaton_factor = value_per_mip/(30000 * 1.6 * 10e-19); // output in mA 
const Double_t rail = 4; // in mA 
const Double_t baseline = 2; // in mA 
const Double_t threshold = 2.5; // in mA

const Double_t noise_mean = 0; 
const Double_t noise_std_dev = 0.01; 



namespace drift{
    const bool cms_approximation = 0; // 
}

namespace frontend{
    const Int_t muxgain_register = 2; // 5 possible values with middle sized register giving gain of 1 mA/mip with +/- 20% 
    const Double_t reference_current = 0.128; // in mA, above which is the gain 
    // const bool peakmode = 0; // 0 for peak mode operation of the APV, 1 for deconvolution mode 
    const Double_t gain = 0.1; // in uA. Need to check which value to use
}
}   // namespace advsnd

#endif   // SHIPLHC_SISENSOR_H_


