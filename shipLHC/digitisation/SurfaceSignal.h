#ifndef SHIPLHC_SURFACESIGNAL_H_
#define SHIPLHC_SURFACESIGNAL_H_

#include "TVector3.h"

#include <iostream>

class SurfaceSignal
{
  public:
    SurfaceSignal()
        : DiffusionArea_(), SurfacePos_(), Amplitude_()
    {
    }

    SurfaceSignal(std::vector<Double_t> DiffusionArea, std::vector<TVector3> SurfacePos, std::vector<Double_t> Amplitude)
        : DiffusionArea_(DiffusionArea)
        , SurfacePos_(SurfacePos)
        , Amplitude_(Amplitude)
    {
    }

    std::vector<Double_t> getDiffusionArea() const { return DiffusionArea_; }
    std::vector<TVector3> getSurfacePos() const { return SurfacePos_; }    
    std::vector<Double_t> getAmplitude() const { return Amplitude_; }    

    void setDiffusionArea(std::vector<Double_t> value) { DiffusionArea_ = value; }
    void setSurfacePos(std::vector<TVector3> value) { SurfacePos_ = value; }
    void setAmplitude(std::vector<Double_t> value) { Amplitude_ = value; }

  private:
    std::vector<Double_t> DiffusionArea_;
    std::vector<TVector3> SurfacePos_;
    std::vector<Double_t> Amplitude_;
};

#endif
