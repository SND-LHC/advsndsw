#ifndef SHIPLHC_CHARGEDIVISION_H_
#define SHIPLHC_CHARGEDIVISION_H_

#include "AdvTargetPoint.h"

#include <iostream>
#include <vector>
// void sumOfTwoNumbers(int a, int b) { std::cout << a << "yo " << b << std::endl; }
class ChargeDivision
{
  public:
    ChargeDivision();
    void ReadPulseShape(std::string PulseFileName);
    std::vector<float> Divide(Int_t detID, const std::vector<AdvTargetPoint*>& V);

  private:
    std::vector<double> PulseValues;
};

#endif   // SHIPLHC_CHARGEDIVISION_H_
