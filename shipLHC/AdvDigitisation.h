#ifndef SHIPLHC_ADVDIGITISATION_H_
#define SHIPLHC_ADVDIGITISATION_H_

#include "AdvTargetPoint.h"
#include "TVector3.h"

#include <iostream>
#include <vector>
using namespace std;

class AdvDigitisation
{
  public:
    AdvDigitisation();
    void digirun(Int_t detID, const std::vector<AdvTargetPoint*>& V);

  protected:
    std::vector<Double_t> EFluct;
    int EFluctSize;
    float segLen;
    std::vector<TVector3> DriftPos;
    std::vector<TVector3> glob_DriftPos;
};
#endif
