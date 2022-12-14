#include "sndMinSciFiConsecutivePlanes.h"

#include "sndSciFiTools.h"

#include "TChain.h"

namespace sndAnalysis{
  minSciFiConsecutivePlanes::minSciFiConsecutivePlanes(TChain * tree) : sciFiBaseCut(tree){
    cutName = "Two or more consecutive SciFi planes";
  }

  bool minSciFiConsecutivePlanes::passCut(){
    initializeEvent();
    
    // For a plane to count, need both planes to have hits
    for (int i = 0; i < hits_per_plane_horizontal.size() - 1; i++){
      if (hits_per_plane_horizontal[i] * hits_per_plane_vertical[i]
	  *hits_per_plane_horizontal[i+1] * hits_per_plane_vertical[i+1] > 0) return true;
    }
    return false;
  }
}	     
