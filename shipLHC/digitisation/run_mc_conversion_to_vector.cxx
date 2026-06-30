#include "AdvPoint.h"
#include "ShipMCTrack.h"

#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"

#include <iostream>
#include <vector>
#include <string>
#include <memory>

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "3 arguments expected but " << argc - 1 << " provided\n";
        std::cerr << "Usage: " << argv[0] << " <input_root_file> <output_root_file> <format: ttree|rntuple>\n";
        return 1;
    }
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);
    std::string output_format(argv[3]);


    auto chain = std::make_unique<TChain>("cbmsim");
    chain->Add(input_path.c_str());
    int nEntries = chain->GetEntries();
    std::cout << nEntries << " entries in TChain\n";

    // Output file and TTree
    auto outfile = TFile::Open(output_path.c_str(), "RECREATE");
    outfile->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD);
    outfile->SetCompressionLevel(5);
    auto outTree = new TTree("cbmsim", "cbmsim");

    // Vectors for output tree
    std::vector<AdvPoint> v_target_points;
    std::vector<ShipMCTrack> v_mc_tracks;

    outTree->Branch("AdvTargetPoint", &v_target_points);
    outTree->Branch("MCTrack", &v_mc_tracks);

    // Input TClonesArrays
    TClonesArray *in_target_points = nullptr;
    TClonesArray *in_mc_tracks = nullptr;

    chain->SetBranchAddress("AdvTargetPoint", &in_target_points);
    chain->SetBranchAddress("MCTrack", &in_mc_tracks);

    for (int i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);

        // Clear vectors and reserve memory
        v_target_points.clear();
        v_target_points.reserve(in_target_points->GetEntries());

        v_mc_tracks.clear();
        v_mc_tracks.reserve(in_mc_tracks->GetEntries());

        // Fill vectors from TClonesArray
        for (int j = 0; j < in_target_points->GetEntries(); ++j) {
            v_target_points.push_back(*(static_cast<AdvPoint*>(in_target_points->At(j))));
        }
        for (int j = 0; j < in_mc_tracks->GetEntries(); ++j) {
            v_mc_tracks.push_back(*(static_cast<ShipMCTrack*>(in_mc_tracks->At(j))));
        }

        outTree->Fill();
    }

    outfile->Write();
    outfile->Close();
}