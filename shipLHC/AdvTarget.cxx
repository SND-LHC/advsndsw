//
//  AdvTarget.cxx
//
//  O. Lantwin and D. Centanni
//  Dec 2023
//

#include "AdvTarget.h"

#include "AdvTargetPoint.h"
#include "FairGeoBuilder.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoNode.h"
#include "FairGeoTransform.h"
#include "FairGeoVolume.h"
#include "FairRootManager.h"
#include "FairVolume.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"
#include "ShipUnit.h"
#include "SiSensor.h"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoGlobalMagField.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoSphere.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"
#include "TGeoUniformMagField.h"
#include "TParticle.h"
#include "TString.h"   // for TString
#include "TVector3.h"
#include "TVirtualMC.h"

#include <ROOT/TSeq.hxx>
#include <iosfwd>     // for ostream
#include <iostream>   // for operator<<, basic_ostream,etc
#include <stddef.h>   // for NULL

using ROOT::TSeq;
using namespace ShipUnit;

AdvTarget::AdvTarget()
    : FairDetector("AdvTarget", "", kTRUE)
    , fTrackID(-1)
    , fVolumeID(-1)
    , fEntryPoint()
    , fMom()
    , fTime(-1.)
    , fLength(-1.)
    , fELoss(-1)
    , fAdvTargetPointCollection(new TClonesArray("AdvTargetPoint"))
{
}

AdvTarget::AdvTarget(const char *name, Bool_t Active, const char *Title)
    : FairDetector(name, true, kAdvSNDTarget)
    , fTrackID(-1)
    , fVolumeID(-1)
    , fEntryPoint()
    , fMom()
    , fTime(-1.)
    , fLength(-1.)
    , fELoss(-1)
    , fAdvTargetPointCollection(new TClonesArray("AdvTargetPoint"))
{
}

AdvTarget::~AdvTarget()
{
    if (fAdvTargetPointCollection) {
        fAdvTargetPointCollection->Delete();
        delete fAdvTargetPointCollection;
    }
}

void AdvTarget::Initialize() { FairDetector::Initialize(); }

// -----   Private method InitMedium
Int_t AdvTarget::InitMedium(const char *name)
{
    static FairGeoLoader *geoLoad = FairGeoLoader::Instance();
    static FairGeoInterface *geoFace = geoLoad->getGeoInterface();
    static FairGeoMedia *media = geoFace->getMedia();
    static FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

    FairGeoMedium *ShipMedium = media->getMedium(name);

    if (!ShipMedium) {
        Fatal("InitMedium", "Material %s not defined in media file.", name);
        return -1111;
    }
    TGeoMedium *medium = gGeoManager->GetMedium(name);
    if (medium != NULL)
        return ShipMedium->getMediumIndex();
    return geoBuild->createMedium(ShipMedium);
}

void AdvTarget::ConstructGeometry()
{
    // Geometry implementation from D. Centanni
    TGeoVolume *top = gGeoManager->GetTopVolume();
    TGeoVolume *detector = gGeoManager->FindVolumeFast("Detector");
    if (!detector)
        LOG(ERROR) << "no Detector volume found ";

    // Materials

    InitMedium("tungstensifon");
    TGeoMedium *tungsten = gGeoManager->GetMedium("tungstensifon");

    InitMedium("Polystyrene");
    TGeoMedium *Polystyrene = gGeoManager->GetMedium("Polystyrene");
    InitMedium("silicon");
    TGeoMedium *Silicon = gGeoManager->GetMedium("silicon");

    Double_t fTargetWallX = conf_floats["AdvTarget/TargetWallX"];
    Double_t fTargetWallY = conf_floats["AdvTarget/TargetWallY"];
    Double_t fTargetWallZ = conf_floats["AdvTarget/TargetWallZ"];
    Double_t fTTZ = conf_floats["AdvTarget/TTZ"];
    Int_t stations = conf_ints["AdvTarget/nTT"];   // Number of TT stations

    TGeoBBox *TargetWall = new TGeoBBox("TargetWall", fTargetWallX / 2., fTargetWallY / 2., fTargetWallZ / 2.);
    TGeoVolume *volTargetWall = new TGeoVolume("volTargetWall", TargetWall, tungsten);
    volTargetWall->SetLineColor(kRed);

    // Silicon tracker module
    //
    // See https://indico.cern.ch/event/1201858/#81-detector-simulation for technical diagrams and renders
    //
    // Passive part
    TGeoBBox *Support = new TGeoBBox("Support", advsnd::module_length / 2, advsnd::module_width / 2, 3.0 * mm / 2);
    TGeoVolume *SupportVolume = new TGeoVolume("SupportVolume", Support, Polystyrene);
    SupportVolume->SetLineColor(kGray);
    // Active part
    TGeoBBox *SensorShape =
        new TGeoBBox("SensorShape", advsnd::sensor_length / 2, advsnd::sensor_width / 2, 0.5 * mm / 2);
    TGeoVolume *SensorVolume = new TGeoVolume("SensorVolumeTarget", SensorShape, Silicon);
    SensorVolume->SetLineColor(kGreen);
    AddSensitiveVolume(SensorVolume);

    // Definition of the target box containing tungsten walls + silicon tracker
    TGeoVolumeAssembly *volAdvTarget = new TGeoVolumeAssembly("volAdvTarget");

    // Positioning calculations, TODO delete once the whole AdvSND apparatus is positioned correctly
    TVector3 EmWall0_survey(5.35 * cm + 42.2 / 2. * cm,
                            17.2 * cm + 42.2 / 2. * cm,
                            288.92 * cm + 10 / 2. * cm);   // Survey position of the centre of the first emulsion wall
    Double_t TargetDiff = 100. * cm - 63.941980 * cm;

    double line_of_sight_offset = (-2.4244059999999976 + 5.203625000000001) * cm;
    detector->AddNode(volAdvTarget,
                      1,
                      new TGeoTranslation(line_of_sight_offset - EmWall0_survey.X() + (40 * cm - 42.2 * cm) / 2.,
                                          EmWall0_survey.Y(),
                                          -TargetDiff + EmWall0_survey.Z() - 60 * cm
                                              - 30 * cm));   // - 60 * cm - 30 * cm to allocate a 150 cm Target

    // For correct detector IDs, the geometry has to be built back to front, from the top-level
    for (auto &&station : TSeq(stations)) {
        TGeoVolumeAssembly *TrackingStation = new TGeoVolumeAssembly("TrackingStation");
        // Each tracking station consists of X and Y planes
        for (auto &&plane : TSeq(advsnd::target::planes)) {
            TGeoVolumeAssembly *TrackerPlane = new TGeoVolumeAssembly("TrackerPlane");
            int i = 0;
            // Each plane consists of 4 modules
            for (auto &&row : TSeq(advsnd::target::rows)) {
                for (auto &&column : TSeq(advsnd::target::columns)) {
                    // Each module in turn consists of two sensors on a support
                    TGeoVolumeAssembly *SensorModule = new TGeoVolumeAssembly("SensorModule");
                    SensorModule->AddNode(SupportVolume, 1);
                    for (auto &&sensor : TSeq(advsnd::sensors)) {
                        int32_t sensor_id =
                            (station << 17) + (plane << 16) + (row << 13) + (column << 11) + (sensor << 10) + 999;
                        SensorModule->AddNode(
                            SensorVolume,
                            sensor_id,
                            new TGeoTranslation(-advsnd::module_length / 2 + 46.95 * mm + advsnd::sensor_length / 2
                                                    + sensor * (advsnd::sensor_length + advsnd::sensor_gap),
                                                0,
                                                +3 * mm / 2 + 0.5 * mm / 2));
                    }
                    TrackerPlane->AddNode(
                        SensorModule,
                        ++i,
                        new TGeoCombiTrans(
                            // Offset modules as needed by row and column
                            TGeoTranslation((column % 2 ? 1 : -1)
                                                * (advsnd::module_length / 2 + advsnd::target::module_column_gap / 2),
                                            (row - 1) * (-advsnd::module_width - advsnd::target::module_row_gap)
                                                - advsnd::target::module_row_gap / 2 + advsnd::module_width / 2,
                                            0),
                            // Rotate every module of the second column
                            TGeoRotation(TString::Format("rot%d", i), 0, 0, column * 180)));
                }
            }
            if (plane == 0) {
                // X-plane
                TrackingStation->AddNode(TrackerPlane, plane, new TGeoTranslation(0, 0, -2 * mm));
            } else if (plane == 1) {
                // Y-plane
                TrackingStation->AddNode(
                    TrackerPlane,
                    plane,
                    new TGeoCombiTrans(TGeoTranslation(0, 0, +2 * mm), TGeoRotation("y_rot", 0, 0, 90)));
            }
        }

        volAdvTarget->AddNode(volTargetWall,
                              station,
                              new TGeoTranslation(0, 0, -fTargetWallZ / 2 + station * fTargetWallZ + station * fTTZ));

        volAdvTarget->AddNode(
            TrackingStation, station, new TGeoTranslation(0, 0, fTTZ / 2 + station * (fTargetWallZ + fTTZ)));
        // new TGeoTranslation(0, 0, (station + 0.5) * fTargetWallZ + (station - 0.5) * 7.5 * mm));
    }
}

Bool_t AdvTarget::ProcessHits(FairVolume *vol)
{
    /** This method is called from the MC stepping */
    // Set parameters at entrance of volume. Reset ELoss.
    if (gMC->IsTrackEntering()) {
        fELoss = 0.;
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fEntryPoint);
        gMC->TrackMomentum(fMom);
    }
    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();

    // Create AdvTargetPoint at exit of active volume
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) {
        if (fELoss == 0.) {
            return kFALSE;
        }

        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();

        TParticle *p = gMC->GetStack()->GetCurrentTrack();
        Int_t pdgCode = p->GetPdgCode();
        TLorentzVector exit_point;
        gMC->TrackPosition(exit_point);
        TLorentzVector Mom;
        gMC->TrackMomentum(Mom);
        Int_t detID = 0;
        gMC->CurrentVolID(detID);
        fVolumeID = detID;
        Double_t xmean = (fEntryPoint.X() + exit_point.X()) / 2.;
        Double_t ymean = (fEntryPoint.Y() + exit_point.Y()) / 2.;
        Double_t zmean = (fEntryPoint.Z() + exit_point.Z()) / 2.;
        AddHit(fTrackID,
               fVolumeID,
               TVector3(xmean, ymean, zmean),
               TVector3(fMom.Px(), fMom.Py(), fMom.Pz()),
               fTime,
               fLength,
               fELoss,
               pdgCode,
               TVector3(exit_point.X(), exit_point.Y(), exit_point.Z()));

        // Increment number of det points in TParticle
        ShipStack *stack = (ShipStack *)gMC->GetStack();
        stack->AddPoint(kAdvSNDTarget);
    }
    return kTRUE;
}

void AdvTarget::GetPosition(Int_t detID, TVector3 &A, TVector3 &B)
{
    // Calculate the detector id as per the geofile, where strips are disrespected
    int strip = (detID) % 1024;                // actual strip ID
    int geofile_detID = detID - strip + 999;   // the det id number needed to read the geometry

    int station = geofile_detID >> 17;
    int plane = (geofile_detID >> 16) % 2;
    int row = (geofile_detID >> 13) % 8;
    int column = (geofile_detID >> 11) % 4;
    int sensor = geofile_detID;
    int sensor_module = advsnd::target::columns * row + 1 + column;

    double global_pos[3];
    double local_pos[3] = {0, 0, 0};
    TString path = TString::Format("/cave_1/"
                                   "Detector_0/"
                                   "volAdvTarget_1/"
                                   "TrackingStation_%d/"
                                   "TrackerPlane_%d/"
                                   "SensorModule_%d/"
                                   "SensorVolumeTarget_%d",
                                   station,
                                   plane,
                                   sensor_module,
                                   sensor);
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    if (nav->CheckPath(path)) {
        nav->cd(path);
    } else {
        LOG(FATAL) << path;
    }
    // Get the corresponding node, which is a sensor made of strips
    TGeoNode *W = nav->GetCurrentNode();
    TGeoBBox *S = dynamic_cast<TGeoBBox *>(W->GetVolume()->GetShape());
    // knowing the strip, get the postion along the sensor
    local_pos[0] = (strip - (advsnd::strips / 2)) * (advsnd::sensor_width / advsnd::strips);
    Double_t top_pos[3] = {local_pos[0], S->GetDY(), 0};
    Double_t bot_pos[3] = {local_pos[0], -(S->GetDY()), 0};
    Double_t global_top_pos[3], global_bot_pos[3];
    nav->LocalToMaster(top_pos, global_top_pos);
    nav->LocalToMaster(bot_pos, global_bot_pos);
    A.SetXYZ(global_top_pos[0], global_top_pos[1], global_top_pos[2]);
    B.SetXYZ(global_bot_pos[0], global_bot_pos[1], global_bot_pos[2]);
}

void AdvTarget::EndOfEvent() { fAdvTargetPointCollection->Clear(); }

void AdvTarget::Register()
{
    /** This will create a branch in the output tree called
        AdvTargetPoint, setting the last parameter to kFALSE means:
        this collection will not be written to the file, it will exist
        only during the simulation.
    */

    FairRootManager::Instance()->Register("AdvTargetPoint", "AdvTarget", fAdvTargetPointCollection, kTRUE);
}
TClonesArray *AdvTarget::GetCollection(Int_t iColl) const
{
    if (iColl == 0) {
        return fAdvTargetPointCollection;
    } else {
        return NULL;
    }
}

void AdvTarget::Reset() { fAdvTargetPointCollection->Clear(); }

AdvTargetPoint *AdvTarget::AddHit(Int_t trackID,
                                  Int_t detID,
                                  TVector3 entrypoint,
                                  TVector3 mom,
                                  Double_t time,
                                  Double_t length,
                                  Double_t eLoss,
                                  Int_t pdgCode,
                                  TVector3 exitpoint)
{
    TClonesArray &clref = *fAdvTargetPointCollection;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) AdvTargetPoint(trackID, detID, entrypoint, mom, time, length, eLoss, pdgCode, exitpoint);
}
