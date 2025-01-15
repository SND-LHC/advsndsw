#ifndef SHIPLHC_ADVSIGNAL_H_
#define SHIPLHC_ADVSIGNAL_H_

#include "TVector3.h"

#include <iostream>

class AdvSignal
{
  public:
    AdvSignal(std::vector<Int_t> Strips, std::vector<Double_t> IntegratedSignal)
        : Strips_(Strips)
        , IntegratedSignal_(IntegratedSignal)
    {
    }

    std::vector<Int_t> getStrips() const { return Strips_; }
    std::vector<Double_t> getIntegratedSignal() const { return IntegratedSignal_; }  

    void setStrips(std::vector<Int_t> value) { Strips_ = value; }
    void setIntegratedSignal(std::vector<Double_t> value) { IntegratedSignal_ = value; }

  private:
    std::vector<Int_t> Strips_;
    std::vector<Double_t> IntegratedSignal_;
};

#endif
