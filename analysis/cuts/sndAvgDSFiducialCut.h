#pragma once

#include "sndMuFilterBaseCut.h"

#include "TChain.h"
#include "MuFilterHit.h"

namespace sndAnalysis {
  class avgDSFiducialCut : public MuFilterBaseCut {
  private :
    double vertical_min, vertical_max, horizontal_min, horizontal_max;
  public :
    avgDSFiducialCut(double vertical_min_cut, double vertical_max_cut, double horizontal_min_cut, double horizontal_max_cut, TChain * tree);
    ~avgDSFiducialCut(){;}

    bool passCut();

  };
};
