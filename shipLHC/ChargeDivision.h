#ifndef SHIPLHC_CHARGEDIVISION_H_
#define SHIPLHC_CHARGEDIVISION_H_

#include "AdvTargetPoint.h"

#include <iostream>
#include <vector>

class ChargeDivision
{
  public:
    ChargeDivision();
    void ReadPulseShape(std::string PulseFileName);
    std::vector<Double_t> Divide(Int_t detID, const std::vector<AdvTargetPoint*>& V);

  private:
    std::vector<double> PulseValues;

  protected:
    Double_t ParticleCharge;
    Double_t ParticleMass;

    float StripPitch = 80e-6;
    float StripWidth = 0.25 * StripPitch;
    Int_t ChargeDivisionsperStrip = 10;
    Int_t NumberofSegments = 0;
    float NumberofStrips = 0;

    float segLen;
};

#endif   // SHIPLHC_CHARGEDIVISION_H_
