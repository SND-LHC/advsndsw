#ifndef SHIPLHC_SURFACESIGNAL_H_
#define SHIPLHC_SURFACESIGNAL_H_

#include "TVector3.h"

#include <iostream>

class SurfaceSignal
{
  public:
    SurfaceSignal(std::vector<Double_t> DiffusionArea, std::vector<TVector3> SurfacePos)
        : DiffusionArea_(DiffusionArea)
        , SurfacePos_(SurfacePos)
    {
    }

    std::vector<Double_t> getDiffusionArea() const { return DiffusionArea_; }
    std::vector<TVector3> getSurfacePos() const { return SurfacePos_; }    

    void setDiffusionArea(std::vector<Double_t> value) { DiffusionArea_ = value; }
    void setSurfacePos(std::vector<TVector3> value) { SurfacePos_ = value; }

  private:
    std::vector<Double_t> DiffusionArea_;
    std::vector<TVector3> SurfacePos_;
};

#endif
