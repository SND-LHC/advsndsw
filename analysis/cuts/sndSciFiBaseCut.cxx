#include "sndSciFiBaseCut.h"

#include "TChain.h"
#include "TClonesArray.h"
#include "sndScifiHit.h"

#include <vector>

namespace sndAnalysis {

TClonesArray* sciFiBaseCut::scifiDigiHitCollection = 0;
TChain* sciFiBaseCut::tree = 0;
unsigned long int sciFiBaseCut::read_entry = -1;

std::vector<int> sciFiBaseCut::hits_per_plane_vertical = std::vector<int>(5, 0);
std::vector<int> sciFiBaseCut::hits_per_plane_horizontal = std::vector<int>(5, 0);

sciFiBaseCut::sciFiBaseCut(TChain* ch)
{
    if (tree == 0) {
        tree = ch;
        scifiDigiHitCollection = new TClonesArray("sndScifiHit", 3000);
        tree->SetBranchAddress("Digi_ScifiHits", &scifiDigiHitCollection);
    }
}

void sciFiBaseCut::initializeEvent()
{
    if (read_entry != tree->GetReadEntry()) {
        read_entry = tree->GetReadEntry();

        // Clear hits per plane vectors
        std::fill(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), 0);
        std::fill(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), 0);

        // Add valid hits to hits per plane vectors
        sndScifiHit* hit;
        TIter hitIterator(scifiDigiHitCollection);

        while ((hit = (sndScifiHit*)hitIterator.Next())) {
            if (hit->isValid()) {
                int sta = hit->GetStation();
                if (hit->isVertical()) {
                    hits_per_plane_vertical[sta - 1]++;
                } else {
                    hits_per_plane_horizontal[sta - 1]++;
                }
            }
        }
    }
}
}   // namespace sndAnalysis
