//
//  AdvTarget.cxx
//
//  S.Ilieva
//  March 2025
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
    TGeoMedium *silicon = gGeoManager->GetMedium("silicon");
    InitMedium("steel");
    TGeoMedium *steel = gGeoManager->GetMedium("steel");

    Double_t fPlateX = conf_floats["AdvTarget/PlateX"];
    Double_t fPlateY = conf_floats["AdvTarget/PlateY"];
    Double_t fPlateZ = conf_floats["AdvTarget/PlateZ"];
    Double_t fPlateFrameX = conf_floats["AdvTarget/PlateFrameX"];
    Double_t fPlateFrameY = conf_floats["AdvTarget/PlateFrameY"];
    Double_t fPlateFrameZ = conf_floats["AdvTarget/PlateFrameZ"];
    Double_t fTTZ = conf_floats["AdvTarget/TTZ"];
    Int_t fNlayers = conf_ints["AdvTarget/nTT"];   // Number of TT layers

    // Definition of the target box containing tungsten walls + silicon tracker
    TGeoVolumeAssembly *volAdvTarget = new TGeoVolumeAssembly("volAdvTarget");

    TGeoBBox *Plate = new TGeoBBox("Plate", fPlateX / 2., fPlateY / 2., fPlateZ / 2. + 1e-10);
    TGeoBBox *PlateAndFrame = new TGeoBBox("PlateAndFrame", fPlateFrameX / 2., fPlateFrameY / 2., fPlateFrameZ / 2.);
    TGeoCompositeShape *PlateFrame = new TGeoCompositeShape("PlateFrame", "PlateAndFrame-Plate");

    TGeoVolume *volPlate = new TGeoVolume("volPlate", Plate, tungsten);
    TGeoVolume *volPlateFrame = new TGeoVolume("volPlateFrame", PlateFrame, steel);
    volPlate->SetLineColor(kGray - 4);
    volPlateFrame->SetLineColor(kGray);
    TGeoVolumeAssembly *volWall = new TGeoVolumeAssembly("volWall");
    volWall->AddNode(volPlate, 1);
    volWall->AddNode(volPlateFrame, 1);

    // Positioning calculations, TODO delete once the whole AdvSND apparatus is positioned correctly
    TVector3 EmWall0_survey(5.35 * cm + 42.2 / 2. * cm,
                            17.2 * cm + 42.2 / 2. * cm,
                            288.92 * cm + 10 / 2. * cm);   // Survey position of the centre of the first emulsion wall
    Double_t TargetDiff = 100. * cm - 63.941980 * cm;

    double line_of_sight_offset = (-2.4244059999999976 + 5.203625000000001) * cm;
    detector->AddNode(volAdvTarget,
                      0,
                      new TGeoTranslation(line_of_sight_offset - EmWall0_survey.X() + (40 * cm - 42.2 * cm) / 2.,
                                          EmWall0_survey.Y(),
                                          -TargetDiff + EmWall0_survey.Z() - 60 * cm
                                              - 30 * cm));   // - 60 * cm - 30 * cm to allocate a 150 cm Target

    // Silicon tracker module
    //
    // See https://indico.cern.ch/event/1452210/contributions/6134180/attachments/2934314/5153496/Abbaneo_Sep2024.pdf
    // for graphical layout
    //
    // Passive part - support with a hole to host the sensors.
    // The support is on 4 legs, which wont be modelled for now. Instead the support_thickness is doubled (2x0.5mm) to
    // approximate.
    TGeoBBox *SupportPlate = new TGeoBBox(
        "SupportPlate", advsnd::module_width / 2, advsnd::module_length / 2, advsnd::support_thickness / 2);
    float support_hole_x_half_dim = (2 * advsnd::sensor_width + advsnd::sensor_gap) / 2.;
    TGeoBBox *SupportHole = new TGeoBBox(
        "SupportHole", support_hole_x_half_dim, advsnd::sensor_length / 2, advsnd::support_thickness / 2 + 1e-10);
    TGeoTranslation *t1 = new TGeoTranslation(
        "t1", -advsnd::module_width / 2 + advsnd::target::module_dead_space_side_small + support_hole_x_half_dim, 0, 0);
    t1->RegisterYourself();
    TGeoCompositeShape *Support = new TGeoCompositeShape("Support", "SupportPlate-SupportHole:t1");

    TGeoVolume *SupportVolume = new TGeoVolume("SupportVolume", Support, Polystyrene);
    SupportVolume->SetLineColor(kGray);
    SupportVolume->SetTransparency(20);
    // Active part  - sensors of Si strips
    // Strips are not included in the geometry since it will become much more computationally heavy
    TGeoBBox *SensorShape =
        new TGeoBBox("SensorShape", advsnd::sensor_width / 2, advsnd::sensor_length / 2, advsnd::sensor_thickness / 2);
    TGeoVolume *SensorVolume = new TGeoVolume("Target_SensorVolume", SensorShape, silicon);
    SensorVolume->SetLineColor(kAzure + 3);
    AddSensitiveVolume(SensorVolume);

    // In detector construction, a station consists of a target(W) plate with modules mounted on the two sides.
    // These modules measure different coordinates e.g. X modules & W & Y modules = station.
    // In the sw the unit is instead the detective layer: modules of the same orientation (almost fully) covering
    // the transverse size of the W plate. A detective layer is formed by the modules in between consecutive target
    // plates. The layer measures X OR Y coordinate! The X and Y layers are alternated in the target volume e.g.
    //           station1         +         station 2            +         station 2           +  ...
    // X modules & W & Y modules  +  Y modules & W & X modules   +  X modules & W & Y modules  +  ...
    // For correct detector IDs, the geometry has to be built back to front, from the bottom up
    for (auto &&layer : TSeq(fNlayers)) {
        // A tracking layer consists of a X OR a Y plane.
        // X plane(vertical strips) is plane 1 and Y plane(horizontal strips) is plane 0.
        // This is done to match the SND-LHC convention, where  plane 1 is vertical
        int plane = layer % advsnd::target::planes;
        TGeoVolumeAssembly *ActiveLayer = new TGeoVolumeAssembly("Target_Layer");
        // Each layer consists of row x columns = 4x2 modules
        int i = 0;
        for (auto &&row : TSeq(advsnd::target::rows)) {
            for (auto &&column : TSeq(advsnd::target::columns)) {
                // Each module in turn consists of two sensors on a support
                TGeoVolumeAssembly *SensorModule = new TGeoVolumeAssembly("SensorModule");
                SensorModule->AddNode(
                    SupportVolume, 1, new TGeoTranslation(0, 0, fTTZ - advsnd::support_thickness / 2));
                for (auto &&sensor : TSeq(advsnd::sensors)) {
                    int32_t sensor_id =
                        (layer << 17) + (plane << 16) + (row << 13) + (column << 11) + (sensor << 10) + 999;
                    SensorModule->AddNode(
                        SensorVolume,
                        sensor_id,
                        new TGeoTranslation(-advsnd::module_width / 2 + advsnd::hcal::module_dead_space_side_small
                                                + advsnd::sensor_width / 2
                                                + sensor * (advsnd::sensor_width + advsnd::sensor_gap),
                                            -advsnd::module_length / 2 + advsnd::hcal::module_dead_space_bottom
                                                + advsnd::sensor_length / 2,
                                            fTTZ - advsnd::support_thickness - advsnd::sensor_thickness / 2));
                }
                // SensorModules one next to the other form the two columns and a little empty space(modules_column_gap)
                // is left in between. When building rows, the overlap btw the upper 4 modules(top tandem) is denoted
                // modules_rows_overlap. The same overlap is applied to the lower 4 modules on the bottom. Note the
                // overlap btw the top tandem and the bottom tandem is different: tandem_modules_rows_overlap.
                ActiveLayer->AddNode(
                    SensorModule,
                    ++i,
                    new TGeoCombiTrans(
                        // Offset modules as needed by row and column
                        TGeoTranslation(
                            (column % 2 ? -1 : 1) * (advsnd::module_width / 2 + advsnd::target::modules_column_gap / 2),
                            (row + 1) * advsnd::module_length / 2
                                + +floor((row + 1) / 2)
                                      * (advsnd::module_length / 2 - advsnd::target::modules_rows_overlap)
                                + floor(row / 2)
                                      * (advsnd::module_length / 2 - advsnd::target::tandem_modules_rows_overlap),
                            fTTZ),   // modules facing one another are offset by the module_thickness, leaving some
                                     // clearance(air) for electronics
                        // Rotate every module of the second column by 180 on z axis and
                        // rotate every other row by 180 deg on y axis to arrive at back-to-front layout
                        TGeoRotation(TString::Format("rot%d", i), 0, row % 2 * 180, column * 180)));
            }   // columns
        }   // rows

        // Alternate target plates followed by a detecting layer.
        // Offsets in Y are adjusted so that the center of the whole detector is positioned at 0
        volAdvTarget->AddNode(volWall, layer, new TGeoTranslation(0, 0, fPlateZ / 2 + layer * (fPlateZ + 2 * fTTZ)));
        if (plane == 0) {
            // Y-plane (horizontal strips along y axis)
            volAdvTarget->AddNode(ActiveLayer,
                                  layer,
                                  new TGeoTranslation(0,
                                                      -fPlateY / 2 + advsnd::target::offset / 2,
                                                      fPlateZ + layer * (fPlateZ + 2 * fTTZ)));
        } else if (plane == 1) {
            // X-plane (vertical strips along x axis)
            volAdvTarget->AddNode(
                ActiveLayer,
                layer,
                new TGeoCombiTrans(TGeoTranslation(2 * advsnd::module_length - advsnd::target::offset,
                                                   -fPlateY / 2 + advsnd::target::offset / 2 + 2 * advsnd::module_length
                                                       - advsnd::target::offset,
                                                   fPlateZ + layer * (fPlateZ + 2 * fTTZ)),
                                   TGeoRotation("y_rot", 0, 0, 90)));
        }
    }   // layers
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

    int layer = geofile_detID >> 17;
    int row = (geofile_detID >> 13) % 8;
    int column = (geofile_detID >> 11) % 4;
    int sensor = geofile_detID;
    int sensor_module = advsnd::target::columns * row + 1 + column;

    double local_pos[3] = {0, 0, 0};
    TString path = TString::Format("/cave_1/"
                                   "Detector_0/"
                                   "volAdvTarget_0/"
                                   "Target_Layer_%d/"
                                   "SensorModule_%d/"
                                   "Target_SensorVolume_%d",
                                   layer,
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
    local_pos[1] = (strip - (advsnd::strips / 2)) * (advsnd::sensor_length / advsnd::strips);
    Double_t left_pos[3] = {S->GetDX(), local_pos[1], 0};
    Double_t right_pos[3] = {-(S->GetDX()), local_pos[1], 0};
    Double_t global_left_pos[3], global_right_pos[3];
    nav->LocalToMaster(left_pos, global_left_pos);
    nav->LocalToMaster(right_pos, global_right_pos);
    A.SetXYZ(global_left_pos[0], global_left_pos[1], global_left_pos[2]);
    B.SetXYZ(global_right_pos[0], global_right_pos[1], global_right_pos[2]);
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
