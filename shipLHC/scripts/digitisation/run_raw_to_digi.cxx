#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include "ROOT/RDataFrame.hxx"

#include "SiStripIO.h"
#include "SiStripRawToDigi.h"

int main(int argc, char* argv[]){
    if (argc != 5) {
        std::cerr << "4 arguments expected but " << argc - 1 << " provided\n";
        std::cerr << "Usage: " << argv[0] << " <input_root_file> <detector_info> <output_root_file> <format: ttree|rntuple>\n";
        return 1;
    }
    std::string input_path(argv[1]);
    std::string detector_info_path(argv[2]);
    std::string output_path(argv[3]);
    std::string output_format(argv[4]);

    // Make format option case insensitve
    std::transform(output_format.begin(), output_format.end(), output_format.begin(), [](unsigned char c) {return std::tolower(c);});

    ROOT::RDF::RSnapshotOptions opts;
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

    auto df = ROOT::RDataFrame("Events", input_path);
    auto df2 = df.Define("FedChannelDigis", SiStripRawToDigi(detector_info_path), {"FEDRawDataCollection_rawDataCollector__LHC."});
    df2.Snapshot("Events", output_path, {"FedChannelDigis"}, opts);
}
