#include "sndMinSciFiHitsCut.h"

#include "sndSciFiTools.h"

#include "TChain.h"

namespace sndAnalysis{
  minSciFiHits::minSciFiHits(int threshold, TChain * tree) : sciFiBaseCut(tree){
    hitThreshold = threshold;
    cutName = "More than "+std::to_string(hitThreshold)+" SciFi hits";
  }

  bool minSciFiHits::passCut(){
    initializeEvent();
    
    if ( getTotalSciFiHits(hits_per_plane_horizontal, hits_per_plane_vertical) < hitThreshold) return false;
    return true;
  }
}	     
