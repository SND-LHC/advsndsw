#ifndef SHIPLHC_FRONTENDDRIVER_H_
#define SHIPLHC_FRONTENDDRIVER_H_

#include "AdvSignal.h"
#include "StripNoise.h"

#include <iostream>
#include <vector>

class FrontendDriver
{
  public:
    FrontendDriver();
    AdvSignal FEDResponse(AdvSignal Signal);
    AdvSignal ADCConversion(AdvSignal ResponseSignal);
    AdvSignal SaturateRange(AdvSignal Signal);
    AdvSignal ZeroSuppressionAlgorithms(AdvSignal Signal);

    void TestingAlgorithm();

  private:
    std::vector<Double_t> ADCcount;

};

#endif  