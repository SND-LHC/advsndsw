#ifndef SHIPLHC_FRONTENDDRIVER_H_
#define SHIPLHC_FRONTENDDRIVER_H_

#include "AdvSignal.h"

#include <iostream>
#include <vector>

class FrontendDriver
{
  public:
    FrontendDriver();
    std::vector<AdvSignal> ADCConversion(std::vector<AdvSignal> ResponseSignal);
    void ZeroSuppressionAlgorithms(AdvSignal Signal);

    void TestingAlgorithm();

  private:
    std::vector<Double_t> ADCcount;
    std::vector<AdvSignal> FEDResponse;

};

#endif  