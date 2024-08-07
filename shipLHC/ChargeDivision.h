#ifndef SHIPLHC_CHARGEDIVISION_H_
#define SHIPLHC_CHARGEDIVISION_H_

#include <iostream>
#include <vector>
// void sumOfTwoNumbers(int a, int b) { std::cout << a << "yo " << b << std::endl; }
class ChargeDivision
{
  public:
    ChargeDivision();
    void ReadPulseShape(std::string test);

  private:
    std::vector<double> PulseValues;
};

#endif   // SHIPLHC_CHARGEDIVISION_H_
