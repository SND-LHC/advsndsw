#ifndef SHIPLHC_FRONTENDDRIVER_H_
#define SHIPLHC_FRONTENDDRIVER_H_

#include "AdvSignal.h"

#include <iostream>
#include <vector>

class FrontedDriver
{
  public:
    FrontedDriver();
    std::vector<AdvSignal> ADCConversion(std::vector<AdvSignal> ResponseSignal);
  private:
    std::vector<Double_t> ADCcount;
    std::vector<AdvSignal> FEDResponse;

};

#endif  