#pragma once

#include <vector>

#include "sndBaseCut.h"

#include "TChain.h"
#include "TClonesArray.h"
#include "sndScifiHit.h"

namespace sndAnalysis {
  
  class sciFiBaseCut : public baseCut {

  private : 
    static TChain * tree;
    static unsigned long int read_entry;

  protected :
    static TClonesArray * scifiDigiHitCollection;

    static std::vector<int> hits_per_plane_vertical;
    static std::vector<int> hits_per_plane_horizontal;

    void initializeEvent();

    sciFiBaseCut(TChain * ch);
    ~sciFiBaseCut(){;}
  };
}
