#include "TPython.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "AdvTargetPoint.h"
#include "AdvTargetHit.h"
#include "AdvMuFilterPoint.h"
#include "AdvMuFilterHit.h"
#include "Hit2MCPoints.h"
//#include "FairLogger.h"
#include "ROOT/RDataFrame.hxx"
#include "digitize.h"

#include <string>
#include <iostream>

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

int main() {
    ROOT::RDF::RSnapshotOptions opts;
    opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kRNTuple;

    auto df = ROOT::RDataFrame("cbmsim", "/afs/cern.ch/work/f/fmei/private/RNTuples-for-advsndsw/data/rntuples-sndLHC.Ntuple-TGeant4.root");

    std::string geometry_path = "/afs/cern.ch/work/f/fmei/private/RNTuples-for-advsndsw/data_multi/geofile_full.Ntuple-TGeant4.root";
    TGeoNavigator* nav = initGeometry(geometry_path);

    auto df2 = df.Define("Digi_AdvTargetHits", advsnd::DigitizePoints<AdvTargetPoint, AdvTargetHit>(nav), {"AdvTargetPoint"});
    auto df3 = df2.Define("Digi_AdvTargetHits2MCPoints", advsnd::LinkPointsToDigi<AdvTargetPoint>(nav), {"AdvTargetPoint"});
    auto df4 = df3.Define("Digi_AdvMuFilterHits", advsnd::DigitizePoints<AdvMuFilterPoint, AdvMuFilterHit>(nav), {"AdvMuFilterPoint"});
    auto df5 = df4.Define("Digi_AdvMuFilterHits2MCPoints", advsnd::LinkPointsToDigi<AdvMuFilterPoint>(nav), {"AdvMuFilterPoint"});
    df5.Snapshot("cbmsim", "/afs/cern.ch/work/f/fmei/private/RNTuples-for-advsndsw/data_multi/rntuples-df-compiled-sndLHC.Ntuple-TGeant4_digCPP.root", {"AdvMuFilterPoint", "AdvTargetPoint", "Digi_AdvTargetHits", "Digi_AdvTargetHits2MCPoints", "Digi_AdvMuFilterHits","Digi_AdvMuFilterHits2MCPoints"}, opts);
}