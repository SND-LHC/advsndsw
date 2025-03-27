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
    InitMedium("silicon");
    TGeoMedium *Silicon = gGeoManager->GetMedium("silicon");
    InitMedium("aluminium");
    TGeoMedium *Aluminum = gGeoManager->GetMedium("aluminium");


    Double_t fTargetWallX = conf_floats["AdvTarget/TargetWallX"];
    Double_t fTargetWallY = conf_floats["AdvTarget/TargetWallY"];
    Double_t fTargetWallZ = conf_floats["AdvTarget/TargetWallZ"];
    Double_t fAluFrameX = conf_floats["AdvTarget/AluFrameX"];
    Double_t fAluFrameY = conf_floats["AdvTarget/AluFrameY"];
    Double_t fAluFrameZ = conf_floats["AdvTarget/AluFrameZ"];
    Double_t fTTX = conf_floats["AdvTarget/TTX"];
    Double_t fTTY = conf_floats["AdvTarget/TTY"];
    Double_t fTTZ = conf_floats["AdvTarget/TTZ"];
    Int_t stations = conf_ints["AdvTarget/nTT"];   // Number of TT stations

    // AdvSND layout August 2024 proposed in https://indico.cern.ch/event/1441955/

    // Definition of the target box containing tungsten walls + silicon tracker (all sensitive plane for now)
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

    TGeoBBox *TargetWall = new TGeoBBox("TargetWall", fTargetWallX/2., fTargetWallY/2., fTargetWallZ/2.);
    TGeoBBox *AluBlock = new TGeoBBox("AluBlock", fAluFrameX/2., fAluFrameY/2., fAluFrameZ/2.-0.001);
    TGeoBBox *Module = new TGeoBBox("Module", fTTX/2., fTTY/2., fTTZ/2.);
    TGeoCompositeShape *AluFrame = new TGeoCompositeShape("AluFrame", "AluBlock-TargetWall");

    TGeoVolume *volTargetWall = new TGeoVolume("volTargetWall", TargetWall, tungsten);
    TGeoVolume *volAluFrame = new TGeoVolume("volAluFrame", AluFrame, Aluminum);
    TGeoVolume *volModule = new TGeoVolume("volModule", Module, Silicon);
    volTargetWall->SetLineColor(kGray-4);
    volAluFrame->SetLineColor(kGray);
    volModule->SetLineColor(kYellow+2);
    AddSensitiveVolume(volModule);
    
    volAdvTarget->AddNode(volModule, 0, new TGeoTranslation(0, 0, fTTZ/2));
    for (auto plane = 1; plane<stations+1; plane++) {
        Double_t station_length = 2*fTTZ+fTargetWallZ;
        Double_t Zshift = fTTZ;
        volAdvTarget->AddNode(volAluFrame, plane, new TGeoTranslation(0, 0, Zshift+fAluFrameZ/2+(plane-1)*station_length));
        volAdvTarget->AddNode(volTargetWall, plane, new TGeoTranslation(0, 0, Zshift+fTargetWallZ/2+(plane-1)*station_length));
        volAdvTarget->AddNode(volModule, plane*10, new TGeoTranslation(0, 0, Zshift+fTargetWallZ+fTTZ/2+(plane-1)*station_length));
        if (plane == stations) continue;
        volAdvTarget->AddNode(volModule, plane*10+1, new TGeoTranslation(0, 0, Zshift+fTargetWallZ+fTTZ+fTTZ/2+(plane-1)*station_length));
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
