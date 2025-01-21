#include "DigiTaskSND.h"   // for Digitization

#include "AdvMuFilterHit.h"
#include "AdvMuFilterPoint.h"
#include "AdvTargetHit.h"
#include "AdvTargetPoint.h"
#include "FairLink.h"            // for FairLink
#include "FairMCEventHeader.h"   // for FairMCEventHeader
#include "FairRootManager.h"     // for FairRootManager
#include "FairRunAna.h"          // for FairRunAna
#include "FairRunSim.h"          // for FairRunSim
#include "Hit2MCPoints.h"        // for linking hits to true MC points
#include "MuFilterHit.h"         // for Muon Filter Hit
#include "MuFilterPoint.h"       // for Muon Filter Point
#include "Scifi.h"               // for SciFi detector
#include "ScifiPoint.h"          // for SciFi Point
#include "SiSensor.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "sndScifiHit.h"   // for SciFi Hit

#include <TClonesArray.h>   // or TClonesArray
#include <TFile.h>
#include <TGenericClassInfo.h>   // for TGenericClassInfo
#include <TMath.h>               // for Sqrt
#include <TROOT.h>
#include <TRandom.h>   // for TRandom, gRandom
#include <algorithm>   // std::sort
#include <iostream>    // for operator<<, basic_ostream, endl
#include <vector>      // std::vector

using namespace std;

DigiTaskSND::DigiTaskSND()
    : FairTask("DigTaskSND")
    , fScifiPointArray(nullptr)
    , fMuFilterPointArray(nullptr)
    , AdvTargetPoints(nullptr)
    , AdvMuFilterPoints(nullptr)
    , fEventHeader(nullptr)
    , fScifiDigiHitArray(nullptr)
    , fMuFilterDigiHitArray(nullptr)
    , fScifiHit2MCPointsArray(nullptr)
    , fMuFilterHit2MCPointsArray(nullptr)
    , AdvTargetHits(nullptr)
    , AdvMuFilterHits(nullptr)
    , AdvTargetHits2MCPoints(nullptr)
    , AdvMuFilterHits2MCPoints(nullptr)
{
}

DigiTaskSND::~DigiTaskSND() {}

InitStatus DigiTaskSND::Init()
{
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman) {
        LOG(FATAL) << "DigiTaskSND::Init: "
                   << "RootManager not instantiated!";
    }

    // Get the SciFi detector and sipm to fibre mapping
    scifi = dynamic_cast<Scifi*>(gROOT->GetListOfGlobals()->FindObject("Scifi"));
    if (scifi) {
        scifi->SiPMmapping();
        fibresSiPM = scifi->GetSiPMmap();
        siPMFibres = scifi->GetFibresMap();
    }

    // Get event header
    // Try classic FairRoot approach first
    fMCEventHeader = static_cast<FairMCEventHeader*>(ioman->GetObject("MCEventHeader."));
    // .. with a safety net for trailing dots mischief
    if (fMCEventHeader == nullptr) {
        fMCEventHeader = static_cast<FairMCEventHeader*>(gROOT->FindObjectAny("MCEventHeader."));
    }
    // Get input MC points
    fScifiPointArray = static_cast<TClonesArray*>(ioman->GetObject("ScifiPoint"));
    fvetoPointArray = static_cast<TClonesArray*>(ioman->GetObject("vetoPoint"));
    fEmulsionPointArray = static_cast<TClonesArray*>(ioman->GetObject("EmulsionDetPoint"));
    fMuFilterPointArray = static_cast<TClonesArray*>(ioman->GetObject("MuFilterPoint"));
    AdvTargetPoints = static_cast<TClonesArray*>(ioman->GetObject("AdvTargetPoint"));
    AdvMuFilterPoints = static_cast<TClonesArray*>(ioman->GetObject("AdvMuFilterPoint"));
    if (!fScifiPointArray and !fMuFilterPointArray) {
        LOG(WARN) << "DigiTaskSND::Init: "
                  << "No Scifi and no MuFilter MC Point array!";
    }
    // copy branches from input file:
    fMCTrackArray = static_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    ioman->Register("MCTrack", "ShipMCTrack", fMCTrackArray, kTRUE);
    if (fvetoPointArray)
        ioman->Register("vetoPoint", "vetoPoints", fvetoPointArray, kTRUE);
    if (fEmulsionPointArray)
        ioman->Register("EmulsionDetPoint", "EmulsionDetPoints", fEmulsionPointArray, kTRUE);
    if (fScifiPointArray)
        ioman->Register("ScifiPoint", "ScifiPoints", fScifiPointArray, kTRUE);
    if (fMuFilterPointArray)
        ioman->Register("MuFilterPoint", "MuFilterPoints", fMuFilterPointArray, kTRUE);
    ioman->Register("AdvTargetPoint", "AdvTargetPoints", AdvTargetPoints, kTRUE);
    ioman->Register("AdvMuFilterPoint", "AdvMuFilterPoints", AdvMuFilterPoints, kTRUE);

    // Event header
    fEventHeader = new SNDLHCEventHeader();
    ioman->Register("EventHeader.", "sndEventHeader", fEventHeader, kTRUE);

    if (fScifiPointArray) {
        // Create and register output array - for SciFi and MuFilter
        fScifiDigiHitArray = new TClonesArray("sndScifiHit");
        ioman->Register("Digi_ScifiHits", "DigiScifiHit_det", fScifiDigiHitArray, kTRUE);
        // Branch containing links to MC truth info
        fScifiHit2MCPointsArray = new TClonesArray("Hit2MCPoints");
        ioman->Register("Digi_ScifiHits2MCPoints", "DigiScifiHits2MCPoints_det", fScifiHit2MCPointsArray, kTRUE);
        fScifiHit2MCPointsArray->BypassStreamer(kTRUE);
    }
    if (fMuFilterPointArray) {
        fMuFilterDigiHitArray = new TClonesArray("MuFilterHit");
        ioman->Register("Digi_MuFilterHits", "DigiMuFilterHit_det", fMuFilterDigiHitArray, kTRUE);
        // Branch containing links to MC truth info
        fMuFilterHit2MCPointsArray = new TClonesArray("Hit2MCPoints");
        ioman->Register(
            "Digi_MuFilterHits2MCPoints", "DigiMuFilterHits2MCPoints_det", fMuFilterHit2MCPointsArray, kTRUE);
        fMuFilterHit2MCPointsArray->BypassStreamer(kTRUE);
    }

    AdvMuFilterHits = new TClonesArray("AdvMuFilterHit");
    ioman->Register("Digi_AdvMuFilterHits", "DigiAdvMuFilterHit_det", AdvMuFilterHits, kTRUE);
    // Branch containing links to MC truth info
    AdvMuFilterHits2MCPoints = new TClonesArray("Hit2MCPoints");
    ioman->Register(
        "Digi_AdvMuFilterHits2MCPoints", "DigiAdvMuFilterHits2MCPoints_det", AdvMuFilterHits2MCPoints, kTRUE);
    AdvMuFilterHits2MCPoints->BypassStreamer(kTRUE);

    AdvTargetHits = new TClonesArray("AdvTargetHit");
    ioman->Register("Digi_AdvTargetHits", "DigiAdvTargetHit_det", AdvTargetHits, kTRUE);
    // Branch containing links to MC truth info
    AdvTargetHits2MCPoints = new TClonesArray("Hit2MCPoints");
    ioman->Register("Digi_AdvTargetHits2MCPoints", "DigiAdvTargetHits2MCPoints_det", AdvTargetHits2MCPoints, kTRUE);
    AdvTargetHits2MCPoints->BypassStreamer(kTRUE);

    return kSUCCESS;
}

void DigiTaskSND::Exec(Option_t* /*opt*/)
{

    if (fScifiDigiHitArray) {
        fScifiDigiHitArray->Clear("C");
        fScifiHit2MCPointsArray->Clear("C");
    }
    if (fMuFilterDigiHitArray) {
        fMuFilterDigiHitArray->Clear("C");
        fMuFilterHit2MCPointsArray->Clear("C");
    }
    AdvTargetHits->Clear("C");
    AdvTargetHits2MCPoints->Clear("C");
    AdvMuFilterHits->Clear("C");
    AdvMuFilterHits2MCPoints->Clear("C");

    // Set event header
    fEventHeader->SetRunId(fMCEventHeader->GetRunID());
    fEventHeader->SetEventNumber(fMCEventHeader->GetEventID());
    fEventHeader->SetBunchType(101);

    // Digitize MC points if any
    if (fMuFilterPointArray)
        digitizeMuFilter();
    if (fScifiPointArray) {
        digitizeScifi();
    }
    digitiseAdvTarget();
    digitiseAdvMuFilter();
}

void DigiTaskSND::digitiseAdvTarget()
{

    int point_index = 0;
    int hit_index = 0;
    std::map<int, std::vector<AdvTargetPoint*>> hit_collector{};
    Hit2MCPoints mc_links;
    std::map<int, std::map<int, double>> mc_points{};
    std::map<int, double> norm{};

    if (!gGeoManager) {
        LOG(FATAL) << "Geofile required to get the position of AdvTargetHits.";
    }
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
    int valid = 0;
    int skip = 0;

    for (auto* ptr : *AdvTargetPoints) {
        auto* point = dynamic_cast<AdvTargetPoint*>(ptr);
        auto detID = point->GetDetectorID();
        if (detID == 0) {
            continue;
        }
        int station = detID / 10;
        int plane = detID % 10;
        auto path = TString::Format("/cave_1/"
                                    "Detector_0/"
                                    "volAdvTarget_1/"
                                    "volModule_%d",
                                    detID
                                    );
        // TODO loop by module?
        if (nav->CheckPath(path)) {
            nav->cd(path);
        } else {
            LOG(FATAL) << path;
        }
        auto* node = nav->GetCurrentNode();
        // Find virtual strip by location
        auto x = point->GetX();
        auto y = point->GetY();
        auto z = point->GetZ();
        double global_pos[3] = {x, y, z};
        double local_pos[3];
        // Move to local coordinates (including rotation) to determine strip
        nav->MasterToLocal(global_pos, local_pos);
        int strip = floor((local_pos[plane] / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2));
        strip = max(0, strip);
        strip = min(advsnd::strips - 1, strip);

        auto detector_id = 10000 * detID + strip;
        // Collect points by virtual strip
        hit_collector[detector_id].emplace_back(point);
        mc_points[detector_id][point_index++] = point->GetEnergyLoss();
        norm[detector_id] += point->GetEnergyLoss();
    }

    for (const auto& [detector_id, points] : hit_collector) {
        // Make one hit per virtual strip (detector ID sensor + strip)
        new ((*AdvTargetHits)[hit_index++]) AdvTargetHit(detector_id, points);
        auto point_map = mc_points[detector_id];
        for (const auto& [point_id, energy_loss] : point_map) {
            mc_links.Add(detector_id, point_id, energy_loss / norm[detector_id]);
        }
    }
    new ((*AdvTargetHits2MCPoints)[0]) Hit2MCPoints(mc_links);
}

void DigiTaskSND::digitiseAdvMuFilter()
{

    int point_index = 0;
    int hit_index = 0;
    std::map<int, std::vector<AdvMuFilterPoint*>> hit_collector{};
    Hit2MCPoints mc_links;
    std::map<int, std::map<int, double>> mc_points{};
    std::map<int, double> norm{};

    if (!gGeoManager) {
        LOG(FATAL) << "Geofile required to get the position of AdvMuFilterHits.";
    }
    TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();

    for (auto* ptr : *AdvMuFilterPoints) {
        auto* point = dynamic_cast<AdvMuFilterPoint*>(ptr);
        auto detID = point->GetDetectorID();
        int station = point->GetStation();
        int plane = point->GetPlane();
        int sensor_module = point->GetModule();
        int sensor = detID;
        auto path = TString::Format("/cave_1/"
                                    "Detector_0/"
                                    "volAdvMuFilter_0/"
                                    "TrackingStation_%d/"
                                    "TrackerPlane_%d/"
                                    "SensorModule_%d/"
                                    "SensorVolumeFilter_%d",
                                    station,
                                    plane,
                                    sensor_module,
                                    sensor);
        // TODO loop by module?
        if (nav->CheckPath(path)) {
            nav->cd(path);
        } else {
            LOG(FATAL) << path;
        }
        auto* node = nav->GetCurrentNode();
        // Find virtual strip by location
        auto x = point->GetX();
        auto y = point->GetY();
        auto z = point->GetZ();
        double global_pos[3] = {x, y, z};
        double local_pos[3];
        // Move to local coordinates (including rotation) to determine strip
        nav->MasterToLocal(global_pos, local_pos);
        int strip = floor((local_pos[0] / (advsnd::sensor_width / advsnd::strips)) + (advsnd::strips / 2));
        strip = max(0, strip);
        strip = min(advsnd::strips - 1, strip);

        auto detector_id = detID - 999 + strip;
        // Collect points by virtual strip
        hit_collector[detector_id].emplace_back(point);
        mc_points[detector_id][point_index++] = point->GetEnergyLoss();
        norm[detector_id] += point->GetEnergyLoss();
    }

    for (const auto& [detector_id, points] : hit_collector) {
        // Make one hit per virtual strip (detector ID sensor + strip)
        new ((*AdvMuFilterHits)[hit_index++]) AdvMuFilterHit(detector_id, points);
        auto point_map = mc_points[detector_id];
        for (const auto& [point_id, energy_loss] : point_map) {
            mc_links.Add(detector_id, point_id, energy_loss / norm[detector_id]);
        }
    }
    new ((*AdvMuFilterHits2MCPoints)[0]) Hit2MCPoints(mc_links);
}

void DigiTaskSND::digitizeScifi()
{
    // a map containing fibreID and vector(list) of points and weights
    map<int, pair<vector<ScifiPoint*>, vector<float>>> hitContainer{};
    Hit2MCPoints mcLinks;
    map<pair<int, int>, double> mcPoints{};
    map<int, double> norm{};
    int globsipmChan{}, detID{};
    int locFibreID{};

    // Fill the map
    for (int k = 0, kEnd = fScifiPointArray->GetEntries(); k < kEnd; k++) {
        ScifiPoint* point = static_cast<ScifiPoint*>(fScifiPointArray->At(k));
        if (!point)
            continue;
        // Collect all hits in same SiPM channel
        detID = point->GetDetectorID();
        locFibreID = detID % 100000;
        // Check if locFibreID in a dead area
        if (siPMFibres.count(locFibreID) == 0)
            continue;
        double dE{};
        float weight{};
        for (auto sipmChan : siPMFibres[locFibreID]) {
            globsipmChan = int(detID / 100000) * 100000 + sipmChan.first;
            weight = sipmChan.second[0];
            hitContainer[globsipmChan].first.push_back(point);
            hitContainer[globsipmChan].second.push_back(weight);
            dE = point->GetEnergyLoss() * weight;
            mcPoints[make_pair(globsipmChan, k)] = dE;
            norm[globsipmChan] += dE;
        }
    }   // End filling map
    int index = 0;
    // Loop over entries of the hitContainer map and collect all hits in same detector element
    for (auto it = hitContainer.begin(); it != hitContainer.end(); it++) {
        new ((*fScifiDigiHitArray)[index])
            sndScifiHit(it->first, hitContainer[it->first].first, hitContainer[it->first].second);
        index++;
        for (auto mcit = mcPoints.begin(); mcit != mcPoints.end(); mcit++) {
            if (it->first == mcit->first.first)
                mcLinks.Add(it->first,
                            mcit->first.second,
                            mcPoints[make_pair(it->first, mcit->first.second)] / norm[it->first]);
        }
    }
    new ((*fScifiHit2MCPointsArray)[0]) Hit2MCPoints(mcLinks);
}

void DigiTaskSND::digitizeMuFilter()
{
    // a map with detID and vector(list) of points
    map<int, vector<MuFilterPoint*>> hitContainer{};
    Hit2MCPoints mcLinks;
    map<pair<int, int>, double> mcPoints{};
    map<int, double> norm{};
    int detID{};

    // Fill the map
    for (int k = 0, kEnd = fMuFilterPointArray->GetEntries(); k < kEnd; k++) {
        MuFilterPoint* point = static_cast<MuFilterPoint*>(fMuFilterPointArray->At(k));
        if (!point)
            continue;
        // Collect all hits in same detector element
        detID = point->GetDetectorID();
        hitContainer[detID].push_back(point);
        mcPoints[make_pair(detID, k)] = point->GetEnergyLoss();
        norm[detID] += point->GetEnergyLoss();
    }
    int index = 0;
    // Loop over entries of the hitContainer map and collect all hits in same detector element
    for (auto it = hitContainer.begin(); it != hitContainer.end(); it++) {
        /*MuFilterHit* aHit = */ new ((*fMuFilterDigiHitArray)[index]) MuFilterHit(it->first, hitContainer[it->first]);
        index++;
        for (auto mcit = mcPoints.begin(); mcit != mcPoints.end(); mcit++) {
            if (it->first == mcit->first.first)
                mcLinks.Add(it->first,
                            mcit->first.second,
                            mcPoints[make_pair(it->first, mcit->first.second)] / norm[it->first]);
        }
    }
    new ((*fMuFilterHit2MCPointsArray)[0]) Hit2MCPoints(mcLinks);
}
