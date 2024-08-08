#ifndef SHIPLHC_ENERGYFLUCTUNIT_H_
#define SHIPLHC_ENERGYFLUCTUNIT_H_

#include "TVector3.h"

#include <iostream>

class EnergyFluctUnit
{
  public:
    EnergyFluctUnit(std::vector<Double_t> Efluct, float segLen, std::vector<TVector3> DriftPos)
        : Efluct_(Efluct)
        , segLen_(segLen)
        , DriftPos_(DriftPos)
    {
    }

    std::vector<Double_t> getEfluct() const { return Efluct_; }
    float getsegLen() const { return segLen_; }
    std::vector<TVector3> getDriftPos() const { return DriftPos_; }

    void setEfluct(std::vector<Double_t> value) { Efluct_ = value; }
    void setsegLen(float value) { segLen_ = value; }
    void setDriftPos(std::vector<TVector3> value) { DriftPos_ = value; }

  private:
    std::vector<Double_t> Efluct_;
    float segLen_;
    std::vector<TVector3> DriftPos_;
};

#endif
