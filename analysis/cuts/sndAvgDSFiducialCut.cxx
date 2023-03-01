#include "sndAvgDSFiducialCut.h"

#include "TChain.h"

namespace sndAnalysis{
  avgDSFiducialCut::avgDSFiducialCut(double vertical_min_cut, double vertical_max_cut, double horizontal_min_cut, double horizontal_max_cut, TChain * tree) : MuFilterBaseCut(tree){
    vertical_min = vertical_min_cut;
    vertical_max = vertical_max_cut;
    horizontal_min = horizontal_min_cut;
    horizontal_max = horizontal_max_cut;
    
    cutName = "Avg DS Ver bar in ["+std::to_string(vertical_min)+","+std::to_string(vertical_max)+"] Hor in ["+std::to_string(horizontal_min)+","+std::to_string(horizontal_max)+"]";

    shortName = "AvgDSbar";
    nbins = std::vector<int>{60, 60};
    range_start = std::vector<double>{0, 60};
    range_end = std::vector<double>{60, 120};
    plot_var = std::vector<double>{-1, -1};
  }

  bool avgDSFiducialCut::passCut(){
    
    double avg_ver = 0.;
    unsigned int n_ver = 0;
    double avg_hor = 0.;
    unsigned int n_hor = 0;

    MuFilterHit* hit;
    TIter hitIterator(muFilterDigiHitCollection);

    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->isValid()){
	if (hit->GetSystem() != 3) continue;

	int x = hit->GetDetectorID() % 1000;

	if (hit->isVertical()){
	  avg_ver += x;
	  n_ver++;
	} else {
	  avg_hor += x;
	  n_hor++;
	}
      }
    }

    if ((n_ver+n_hor) == 0) {
      plot_var[0] = -1;
      plot_var[1] = -1;
      return false;
    }
    
    if (n_ver) {
      avg_ver /= n_ver;
      plot_var[0] = avg_ver;
    } else {
      plot_var[0] = -1;
    }

    if (n_hor) {
      avg_hor /= n_hor;
      plot_var[1] = avg_hor;
    } else {
      plot_var[1] = -1;
    }

    if (n_ver == 0) return false;
    if (n_hor == 0) return false;

    if (avg_hor < horizontal_min) return false;
    if (avg_hor > horizontal_max) return false;
    if (avg_ver < vertical_min) return false;
    if (avg_ver > vertical_max) return false;

    return true;
  }
}	     
