#include "sndSciFiStationCut.h"

#include "sndSciFiTools.h"

#include "TChain.h"

namespace sndAnalysis{
  sciFiStationCut::sciFiStationCut(float threshold, std::vector<int> excluded_stations, TChain * tree) : sciFiBaseCut(tree){
    fractionThreshold = threshold;
    stations_to_exclude = std::vector(excluded_stations);
    cutName = "Exclude stations";
    for (int sta : excluded_stations){
      cutName += " "+std::to_string(sta);
    }
    cutName += ". Threshold "+std::to_string(fractionThreshold);
  }

  bool sciFiStationCut::passCut(){
    initializeEvent();
    
    int station = findStation(hits_per_plane_horizontal, hits_per_plane_vertical, fractionThreshold);
    
    if (std::find(stations_to_exclude.begin(), stations_to_exclude.end(), station) == stations_to_exclude.end()){
      return true;
    } else {
      return false;
    }
  }
}	     
