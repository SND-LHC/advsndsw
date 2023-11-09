#include "sndUSQDCCut.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"
#include "TString.h"

#include <vector>
#include <map>
#include <numeric>

namespace sndAnalysis {

  USQDCCut::USQDCCut(float threshold, TChain * ch) : MuFilterBaseCut(ch) {
    qdc_threshold = threshold;
    cutName = "Total US QDC > "+std::to_string(qdc_threshold);

    shortName = "USQDC";
    nbins = std::vector<int>{100};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{10000};
    plot_var = std::vector<double>{-1};
  }

  bool USQDCCut::passCut(){
    MuFilterHit * hit;
    TIter hitIterator(muFilterDigiHitCollection);

    float totQDC = 0.;

    bool ds = false;
    std::vector<bool> us = std::vector<bool>(5, false); 
    
    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->GetSystem() == 2) {
	for (const auto& [key, value] : hit->GetAllSignals()) {
	  totQDC += value;
	}
      }
    }
    plot_var[0] = totQDC;
    if (totQDC >= qdc_threshold) return true;
    return false;
  }
}
