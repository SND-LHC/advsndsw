#ifndef SHIPLHC_SIDIGIPARAMETERS_H_
#define SHIPLHC_SIDIGIPARAMETERS_H_

/* Header-only to implement common code for CMS TOB */
#include "ShipUnit.h"

namespace stripsensor {

const bool peakmode = 0;
const Double_t sampling = 25; 
const std::string APVpeakpulse = "advsndsw/shipLHC/data/APVShapePeak_default.txt";
const std::string APVdecopulse = "advsndsw/shipLHC/data/APVShapeDeco_default.txt";

const Double_t value_per_mip = 1; // in mA
const Double_t amplificaton_factor = value_per_mip/(30000 * 1.6 * 10e-19); // output in mA 
const Double_t rail = 4; // in mA 
const Double_t baseline = 2; // in mA 
const Double_t threshold = 2.5; // in mA
const Double_t mode = 20 ; // in MHz - the frequency that output is cycled 

const Double_t noise_mean = 0; 
const Double_t noise_std_dev = 0.01; 

namespace chargedivision{
    const Int_t ChargeDivisionsperStrip = 10;
    const float StripPitch = 120e-4;
}



namespace drift{
    const bool cms_approximation = 0; // 
    const Double_t module_thickness = 0.05; //check value
    const int depletion_voltage = 170; 
    const int applied_voltage = 300; 
    const int charge_mobility = 480;
    const Double_t chargedistributionRMS = 6.5e-10;
    const int temperature = 300; 
    const Double_t perGeV = 3.61; 

}

namespace inducedcharge{
    const Int_t NSigma = 3; 
    const Double_t strip_width = 30e-4; //in cm  
    const Double_t strip_pitch = 120e-4; //in cm 
}

namespace frontend{
    const Int_t muxgain_register = 2; // 5 possible values with middle sized register giving gain of 1 mA/mip with +/- 20% 
    const Double_t reference_current = 0.128; // in mA, above which is the gain 
    // const bool peakmode = 0; // 0 for peak mode operation of the APV, 1 for deconvolution mode 
    const Double_t gain = 0.1; // in uA. Need to check which value to use
    const bool write_digi_to_text = 1; 
    const Int_t ElectronperADC = 250; 
    const Int_t StripNoise = 5; 
}
}   // namespace advsnd

#endif   // SHIPLHC_SISENSOR_H_


