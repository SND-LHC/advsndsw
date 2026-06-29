#include "TPython.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "AdvPoint.h"
#include "AdvHit.h"
#include "Hit2MCPoints.h"
#include "ROOT/RDataFrame.hxx"
#include "run_mc_digitisation.h"

#include <string>
#include <iostream>
#include <stdlib.h>

TGeoNavigator* initGeometry(const std::string& geometry_path)
{
    TPython::Exec("import SndlhcGeo");
    TPython::Exec(("SndlhcGeo.GeoInterface('" + geometry_path + "')").c_str());

    if (!gGeoManager) {
        std::cerr << "Geofile required\n";
    }
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
    if (!nav) {
        std::cerr << "Failed to get TGeoNavigator from gGeoManager.\n";
    }
    return nav;
}

int main(int argc, char *argv[]) {
        if (argc != 5) {
        std::cerr << "4 arguments expected but " << argc - 1 << " provided\n";
        std::cerr << "Usage: raw_to_digi <input_root_file> <output_root_file> <geometry_file> <format: ttree|rntuple>\n";
        return 1;
    }
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);
    std::string geo_path(argv[3]);
    std::string output_format(argv[4]);

    // Make format option case insensitve
    std::transform(output_format.begin(), output_format.end(), output_format.begin(), [](unsigned char c) {return std::tolower(c);});

    ROOT::RDF::RSnapshotOptions opts;
    opts.fVector2RVec = false;
    opts.fCompressionAlgorithm = ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD;
    opts.fCompressionLevel = 5;
    if (output_format == "rntuple") {
        opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kRNTuple;
        std::cout << "Using RNTuple output format\n";
    }
    else if (output_format == "ttree") {
        opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kTTree;
        std::cout << "Using TTree output format\n";
    }
    else {
        std::cerr << "Unknown output format: " << output_format << "\n";
        std::cerr << "Supported formats: ttree, rntuple\n";
        return 1;
    }

    TGeoNavigator* nav = initGeometry(geo_path);

    auto df = ROOT::RDataFrame("cbmsim", input_path);
    std::cout << df.GetColumnType("AdvPoint") << std::endl;

    auto df2 = df.Define("Digi_AdvHits", advsnd::DigitizePoints(nav), {"AdvPoint"})
        .Define("Digi_AdvHits2MCPoints", advsnd::LinkPointsToDigi(nav), {"AdvPoint"});
    df2.Snapshot("cbmsim", output_path, {"MCTrack", "AdvPoint", "Digi_AdvHits", "Digi_AdvHits2MCPoints"}, opts);
}