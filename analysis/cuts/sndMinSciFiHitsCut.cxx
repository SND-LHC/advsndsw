#include "sndMinSciFiHitsCut.h"

#include "sndSciFiTools.h"

#include "TChain.h"

namespace sndAnalysis{
  minSciFiHits::minSciFiHits(int threshold, TChain * tree) : sciFiBaseCut(tree){
    hitThreshold = threshold;
    cutName = "More than "+std::to_string(hitThreshold)+" SciFi hits";
    shortName = "SciFiMinHits";
    nbins = std::vector<int>{1536};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{1536};
    plot_var = std::vector<double>{-1};
  }

  bool minSciFiHits::passCut(){
    initializeEvent();
    plot_var[0] = getTotalSciFiHits(hits_per_plane_horizontal, hits_per_plane_vertical);
    if ( plot_var[0] < hitThreshold) return false;
    return true;
  }
}	     
