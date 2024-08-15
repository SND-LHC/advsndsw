#ifndef SHIPLHC_ENERGYFLUCTUNIT_H_
#define SHIPLHC_ENERGYFLUCTUNIT_H_

#include "TVector3.h"

#include <iostream>

class EnergyFluctUnit
{
  public:
    EnergyFluctUnit(std::vector<Double_t> Efluct, float segLen, std::vector<TVector3> DriftPos, std::vector<TVector3> glob_DriftPos)
        : Efluct_(Efluct)
        , segLen_(segLen)
        , DriftPos_(DriftPos)
        , glob_DriftPos_(glob_DriftPos)
    {
    }

    std::vector<Double_t> getEfluct() const { return Efluct_; }
    float getsegLen() const { return segLen_; }
    std::vector<TVector3> getDriftPos() const { return DriftPos_; }
    std::vector<TVector3> getglobDriftPos() const { return glob_DriftPos_; }
    

    void setEfluct(std::vector<Double_t> value) { Efluct_ = value; }
    void setsegLen(float value) { segLen_ = value; }
    void setDriftPos(std::vector<TVector3> value) { DriftPos_ = value; }
    void setglobDriftPos(std::vector<TVector3> value) { glob_DriftPos_ = value; }

  private:
    std::vector<Double_t> Efluct_;
    float segLen_;
    std::vector<TVector3> DriftPos_;
    std::vector<TVector3> glob_DriftPos_;
};

#endif
