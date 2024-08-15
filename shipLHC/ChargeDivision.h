#ifndef SHIPLHC_CHARGEDIVISION_H_
#define SHIPLHC_CHARGEDIVISION_H_

#include "AdvTargetPoint.h"
#include "EnergyFluctUnit.h"
#include "TVector3.h"

#include <iostream>
#include <vector>

class ChargeDivision
{
  public:
    ChargeDivision();
    void ReadPulseShape(std::string PulseFileName);
    EnergyFluctUnit Divide(Int_t detID, const std::vector<AdvTargetPoint*>& V);
    TVector3 DriftDir(TVector3 EntryPoint, TVector3 ExitPoint, float length);
    TVector3 getLocal(Int_t detID,TVector3 point);

  private:
    std::vector<double> PulseValues;

  protected:
    Double_t ParticleCharge;
    Double_t ParticleMass;

    float StripPitch = 120e-4;
    //float StripWidth = 0.25 * StripPitch;
    Int_t ChargeDivisionsperStrip = 10;
    Int_t NumberofSegments = 0;
    float NumberofStrips = 0;

    float segLen;
};

#endif   // SHIPLHC_CHARGEDIVISION_H_
