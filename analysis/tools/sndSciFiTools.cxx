#include "sndSciFiTools.h"

#include <numeric>
#include <algorithm>

#include "TClonesArray.h"
#include "sndScifiHit.h"

void sndAnalysis::getSciFiHitsPerStation(TClonesArray * digiHits, std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits){

  // Clear hits per plane vectors
  std::fill(horizontal_hits.begin(), horizontal_hits.end(), 0);
  std::fill(vertical_hits.begin(), vertical_hits.end(), 0);

  // Add valid hits to hits per plane vectors
  sndScifiHit * hit;
  TIter hitIterator(digiHits);

  while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
    if (hit->isValid()){
      int sta = hit->GetStation();
      if (hit->isVertical()){
	vertical_hits[sta-1]++;
      } else {
	horizontal_hits[sta-1]++;
      }
    }
  }
}

int sndAnalysis::getTotalSciFiHits(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits){
  return std::accumulate(horizontal_hits.begin(),
			 horizontal_hits.end(),
			 std::accumulate(vertical_hits.begin(),
					 vertical_hits.end(), 0));
}

int sndAnalysis::getTotalSciFiHits(TClonesArray * digiHits){
  std::vector<int> horizontal_hits = std::vector<int>(5);
  std::vector<int> vertical_hits = std::vector<int>(5);
    
  getSciFiHitsPerStation(digiHits, horizontal_hits, vertical_hits);

  return getTotalSciFiHits(horizontal_hits, vertical_hits);
}

std::vector<float> sndAnalysis::getFractionalHitsPerPlane(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits){
    
  int total_hits = getTotalSciFiHits(horizontal_hits, vertical_hits);
    
  std::vector<float> fractional_hits_per_station = std::vector<float>(horizontal_hits.size());
    
  std::transform(horizontal_hits.begin(), horizontal_hits.end(),
		 vertical_hits.begin(),
		 fractional_hits_per_station.begin(),
		 [&total_hits](const auto& hor, const auto& ver){
		   return ((float) hor + ver)/total_hits;
		 });
    
  return fractional_hits_per_station;
    
}

std::vector<float> sndAnalysis::getFractionalHitsPerPlane(TClonesArray * digiHits){
  std::vector<int> horizontal_hits = std::vector<int>(5);
  std::vector<int> vertical_hits = std::vector<int>(5);
    
  getSciFiHitsPerStation(digiHits, horizontal_hits, vertical_hits);
    
  return getFractionalHitsPerPlane(horizontal_hits, vertical_hits);
}
  
  
int sndAnalysis::findStation(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits, float threshold){

  std::vector<float> frac = getFractionalHitsPerPlane(horizontal_hits, vertical_hits);
    
  std::vector<float> frac_sum = std::vector<float>(frac.size());
    
  std::partial_sum(frac.begin(), frac.end(), frac_sum.begin());

  std::vector<float>::iterator station = std::find_if(frac_sum.begin(), frac_sum.end(),
						      [&threshold](const auto& f) {
							return f > threshold;
						      });
    
  return station - frac_sum.begin() + 1;
    
}

int sndAnalysis::findStation(TClonesArray * digiHits, float threshold){
  std::vector<int> horizontal_hits = std::vector<int>(5);
  std::vector<int> vertical_hits = std::vector<int>(5);
    
  getSciFiHitsPerStation(digiHits, horizontal_hits, vertical_hits);

  return findStation(horizontal_hits, vertical_hits, threshold);
}

