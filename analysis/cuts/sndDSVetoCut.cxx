#include "sndDSVetoCut.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"

#include <vector>
#include <numeric>

namespace sndAnalysis {

  DSVetoCut::DSVetoCut(TChain * ch) : MuFilterBaseCut(ch) {
    cutName = "Remove events with hits in the last (hor) and two last (ver) DS planes";
    
    shortName = "DSVetoCut";
    nbins = std::vector<int>{180};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{180};
    plot_var = std::vector<double>{-1};

  }

  bool DSVetoCut::passCut(){
    MuFilterHit * hit;
    TIter hitIterator(muFilterDigiHitCollection);

    //    bool ds = false;
    //    std::vector<bool> us = std::vector<bool>(5, false); 

    double n_hits = 0;
    
    bool ret = true;
    
    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->GetSystem() == 3) { // DS
	if (hit->GetPlane() >= 3) {
	  ret = false; 
	  n_hits+=1;
	}
      }
    }
    
    plot_var[0] = n_hits;
    return ret;
  }
}
