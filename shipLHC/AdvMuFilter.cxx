//
//  AdvMuFilter.cxx
//
//  S.Ilieva
//  March 2025
//
//  by D. Centanni and O. Lantwin
//  2024
//

#include "AdvMuFilter.h"

#include "AdvMuFilterPoint.h"
#include "FairGeoBuilder.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoNode.h"
#include "FairGeoTransform.h"
#include "FairGeoVolume.h"
#include "FairRootManager.h"
#include "FairRun.h"   // for FairRun
#include "FairRun.h"
#include "FairRuntimeDb.h"   // for FairRuntimeDb
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"
#include "ShipUnit.h"
#include "SiSensor.h"
#include "TClonesArray.h"
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
#include "TList.h"       // for TListIter, TList (ptr only)
#include "TObjArray.h"   // for TObjArray
#include "TParticle.h"
#include "TString.h"   // for TString
#include "TVector3.h"
#include "TVirtualMC.h"

#include <ROOT/TSeq.hxx>
#include <cstring>
#include <iosfwd>     // for ostream
#include <iostream>   // for operator<<, basic_ostream,etc
#include <stddef.h>   // for NULL
#include <string.h>

using ROOT::TSeq;
using std::cout;
using std::endl;
using std::string;
using std::to_string;

using namespace ShipUnit;

AdvMuFilter::AdvMuFilter()
    : FairDetector("AdvMuFilter", "", kTRUE)
    , fTrackID(-1)
    , fVolumeID(-1)
    , fPos()
    , fMom()
    , fTime(-1.)
    , fLength(-1.)
    , fELoss(-1)
    , fAdvMuFilterPointCollection(new TClonesArray("AdvMuFilterPoint"))
{
}

AdvMuFilter::AdvMuFilter(const char *name, Bool_t Active, const char *Title)
    : FairDetector(name, true, kAdvSNDMuFilter)
    , fTrackID(-1)
    , fVolumeID(-1)
    , fPos()
    , fMom()
    , fTime(-1.)
    , fLength(-1.)
    , fELoss(-1)
    , fAdvMuFilterPointCollection(new TClonesArray("AdvMuFilterPoint"))
{
}

AdvMuFilter::~AdvMuFilter()
{
    if (fAdvMuFilterPointCollection) {
        fAdvMuFilterPointCollection->Delete();
        delete fAdvMuFilterPointCollection;
    }
}

void AdvMuFilter::Initialize() { FairDetector::Initialize(); }

// -----   Private method InitMedium
Int_t AdvMuFilter::InitMedium(const char *name)
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

void AdvMuFilter::ConstructGeometry()
{
    // Geometry implementation from D. Centanni
    TGeoVolume *top = gGeoManager->GetTopVolume();
    TGeoVolume *detector = gGeoManager->FindVolumeFast("Detector");
    if (!detector)
        LOG(ERROR) << "no Detector volume found ";

    // Materials

    InitMedium("iron");
    TGeoMedium *Fe = gGeoManager->GetMedium("iron");

    InitMedium("CoilCopper");
    TGeoMedium *Cu = gGeoManager->GetMedium("CoilCopper");

    InitMedium("Polystyrene");
    TGeoMedium *Polystyrene = gGeoManager->GetMedium("Polystyrene");

    InitMedium("aluminium");
    TGeoMedium *Al = gGeoManager->GetMedium("aluminium");

    InitMedium("silicon");
    TGeoMedium *silicon = gGeoManager->GetMedium("silicon");

    InitMedium("steel");
    TGeoMedium *steel = gGeoManager->GetMedium("steel");

    Double_t fMuonSysPlaneX = conf_floats["AdvMuFilter/MuonSysPlaneX"];
    Double_t fMuonSysPlaneY = conf_floats["AdvMuFilter/MuonSysPlaneY"];
    Double_t fMuonSysPlaneZ = conf_floats["AdvMuFilter/MuonSysPlaneZ"];
    Double_t fFeX = conf_floats["AdvMuFilter/FeX"];
    Double_t fFeY = conf_floats["AdvMuFilter/FeY"];
    Double_t fFeZ = conf_floats["AdvMuFilter/FeZ"];
    Double_t fFeGap = conf_floats["AdvMuFilter/FeGap"];
    Double_t fCurvRadius = conf_floats["AdvMuFilter/CurvRadius"];
    Int_t fNlayers = advsnd::hcal::n_XY_layers + advsnd::hcal::n_X_layers;
    Double_t fCoilZ = conf_floats["AdvMuFilter/CoilZ"];
    Double_t fCoilY = conf_floats["AdvMuFilter/CoilY"];
    Double_t fCoilX = fMuonSysPlaneX + 2 * fFeGap - 2 * (fCoilZ + fCurvRadius);
    Double_t fFeYokeX = (fFeX - fMuonSysPlaneX - 2 * fFeGap) / 2.;
    Double_t fFeYokeY = (fFeY - fMuonSysPlaneY) / 2.;

    LOG(INFO) << " CoilX: " << fCoilX << " cm" << endl;
    LOG(INFO) << " FeYokeX: " << fFeYokeX << " cm" << endl;
    LOG(INFO) << " FeYokeY: " << fFeYokeY << " cm" << endl;

    TGeoVolumeAssembly *volAdvMuFilter = new TGeoVolumeAssembly("volAdvMuFilter");

    // Moving the detector forward in the Z direction adding:
    //	-	354.362 which is the last target detector plane midpoint position
    // 	-	3 is the target detector plane half Z-dimension
    //	-	1 is the clearance between volAdvTarget and volAdvMuilter

    Double_t fTargetWallX = 50.;   // cm
    TVector3 EmWall0_survey(
        5.35 + 42.2 / 2. - (fTargetWallX - 42.2) / 2., 17.2 + 42.2 / 2., 288.92 + 10 / 2. + 100);   // cm

    detector->AddNode(volAdvMuFilter,
                      0,
                      new TGeoTranslation(-2.4244059999999976 - EmWall0_survey.X() - 1.1,
                                          38.3,
                                          354.862 + 3 + 1 + fFeZ / 2. - 41.895793 + 1. - 3.854227000000008
                                              + 3.7497190000000046
                                              - 64.44596200000001));   // hardcoded, try to find and elegant solution

    // Iron and Detector's shapes
    TGeoBBox *FeBlock = new TGeoBBox("FeBlock", fFeX / 2., fFeY / 2., fFeZ / 2.);
    TGeoBBox *FeGap = new TGeoBBox("FeGap", (2 * fFeGap + fMuonSysPlaneX) / 2, fMuonSysPlaneY / 2., fFeZ / 2. + 1e-10);
    TGeoCompositeShape *FeSlab = new TGeoCompositeShape("FeSlab", "FeBlock-FeGap");
    TGeoBBox *MagnetizedFe = new TGeoBBox("MagnetizedFe", fMuonSysPlaneX / 2., fMuonSysPlaneY / 2., fFeZ / 2.);

    // Coil shapes
    TGeoBBox *FrontCoil = new TGeoBBox("FrontCoil", fCoilX / 2., fCoilY / 2., fCoilZ / 2.);
    TGeoBBox *LatCoil = new TGeoBBox("LatCoil", fCoilZ / 2., fCoilY / 2., fNlayers * (fFeZ + 2 * fMuonSysPlaneZ) / 2.);
    TGeoTubeSeg *CurvCoil = new TGeoTubeSeg("CurvCoil", fCurvRadius, fCurvRadius + fCoilZ, fCoilY / 2., 0, 90);

    TGeoCombiTrans *Left =
        new TGeoCombiTrans(TGeoTranslation(-fCoilX / 2., 0, -(fCoilZ + fCurvRadius) / 2. - fCurvRadius / 2.),
                           TGeoRotation("rot1", 0, 90, 90));
    Left->SetName("Left");
    Left->RegisterYourself();
    TGeoCombiTrans *Right =
        new TGeoCombiTrans(TGeoTranslation(+fCoilX / 2., 0, -(fCoilZ + fCurvRadius) / 2. - fCurvRadius / 2.),
                           TGeoRotation("rot2", 0, 90, 0));
    Right->SetName("Right");
    Right->RegisterYourself();
    TGeoCompositeShape *FrontCoilShape =
        new TGeoCompositeShape("FrontCoilShape", "FrontCoil+(CurvCoil:Left)+(CurvCoil:Right)");

    TGeoVolume *volFrontCoil = new TGeoVolume("volFrontCoil", FrontCoilShape, Cu);
    volFrontCoil->SetLineColor(kOrange + 1);
    volAdvMuFilter->AddNode(
        volFrontCoil, 0, new TGeoCombiTrans(TGeoTranslation(0, 0, fCoilZ / 2.), TGeoRotation("rot3", 0, 180, 0)));
    volAdvMuFilter->AddNode(
        volFrontCoil,
        1,
        new TGeoTranslation(0, 0, fCoilZ + 2 * fCurvRadius + fCoilZ / 2. + fNlayers * (fFeZ + 2 * fMuonSysPlaneZ)));
    TGeoTranslation *LatLeft = new TGeoTranslation(-fCoilX / 2. - fCoilZ / 2. - fCurvRadius, 0, 0);
    LatLeft->SetName("LatLeft");
    LatLeft->RegisterYourself();
    TGeoTranslation *LatRight = new TGeoTranslation(+fCoilX / 2. + fCoilZ / 2. + fCurvRadius, 0, 0);
    LatRight->SetName("LatRight");
    LatRight->RegisterYourself();
    TGeoCompositeShape *LatCoilShape = new TGeoCompositeShape("LatCoilShape", "(LatCoil:LatLeft)+(LatCoil:LatRight)");
    TGeoVolume *volLatCoil = new TGeoVolume("volLatCoil", LatCoilShape, Cu);
    volLatCoil->SetLineColor(kOrange + 1);
    volAdvMuFilter->AddNode(
        volLatCoil, 0, new TGeoTranslation(0, 0, fCoilZ + fCurvRadius + (fNlayers / 2.) * (fFeZ + 2 * fMuonSysPlaneZ)));

    TGeoVolume *volFeSlab = new TGeoVolume("volFeSlab", FeSlab, Fe);
    TGeoVolume *volMagnetizedFe = new TGeoVolume("volMagnetizedFe", MagnetizedFe, Fe);

    Double_t fField = conf_floats["AdvMuFilter/Field"];
    LOG(INFO) << " Magnetic field: " << fField / 10. << " Tesla" << endl;
    TGeoUniformMagField *magField = new TGeoUniformMagField(0, fField, 0);
    TGeoGlobalMagField::Instance()->SetField(magField);
    volMagnetizedFe->SetField(magField);

    volFeSlab->SetLineColor(kGreen - 4);
    volMagnetizedFe->SetLineColor(kGreen);

    // No excavation case with smaller HCAL layers. The silicon strip modules are identical to the target ones.
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
    TGeoVolume *SensorVolume = new TGeoVolume("HCAL_SensorVolume", SensorShape, silicon);
    SensorVolume->SetLineColor(kAzure + 3);
    AddSensitiveVolume(SensorVolume);

    // Useful dimensions when placing detector elements
    float Zshift = fMuonSysPlaneZ + fCoilZ + fCurvRadius;
    float layer_DX = (2 * advsnd::module_width + advsnd::hcal::modules_column_gap) / 2;
    float layer_DY = (4 * advsnd::module_length - 3 * advsnd::hcal::modules_rows_overlap) / 2;
    float Fe_and_layer_z = fFeZ + 2 * fMuonSysPlaneZ;

    for (auto &&layer : TSeq(fNlayers)) {
        // if (i == fNlayers - 1)
        //     continue;

        TGeoVolumeAssembly *ActiveLayer = new TGeoVolumeAssembly("HCAL_Layer");
        int plane = (layer < advsnd::hcal::n_XY_layers) ? layer % advsnd::hcal::planes : 1;
        // Each plane consists of row x columns = 4x2 modules
        int i = 0;
        for (auto &&row : TSeq(advsnd::hcal::rows)) {
            for (auto &&column : TSeq(advsnd::hcal::columns)) {
                // Each module in turn consists of two sensors on a support
                TGeoVolumeAssembly *SensorModule = new TGeoVolumeAssembly("SensorModule");
                SensorModule->AddNode(
                    SupportVolume, 1, new TGeoTranslation(0, 0, fMuonSysPlaneZ - advsnd::support_thickness / 2));
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
                                            fMuonSysPlaneZ - advsnd::support_thickness - advsnd::sensor_thickness / 2));
                }
                ActiveLayer->AddNode(
                    SensorModule,
                    ++i,
                    new TGeoCombiTrans(
                        // Offset modules as needed by column and rotate
                        TGeoTranslation(
                            (column % 2 ? -1 : 1) * (advsnd::module_width / 2 + advsnd::hcal::modules_column_gap / 2),
                            (row + 1) * advsnd::module_length / 2
                                + +row * (advsnd::module_length / 2 - advsnd::hcal::modules_rows_overlap),
                            fMuonSysPlaneZ),   // modules facing one another are offset by the module_thickness, leaving
                                               // some clearance(air) for electronics
                        // Rotate every module of the first column by 180 on z axis and
                        // rotate every other row by 180 deg on y axis to arrive at back-to-front layout
                        TGeoRotation(TString::Format("rot%d", i), 0, (row % 2 ? 180 : 0), (column % 2 ? 180 : 0))));
            }   // columns
        }   // rows
        // Alternate iron slabs followed by a detecting layer
        // Offsets in Y are adjusted so that the center of the whole detector is positioned at 0
        volAdvMuFilter->AddNode(
            volFeSlab, layer, new TGeoTranslation(0, 0, Zshift + fFeZ / 2 + layer * Fe_and_layer_z));
        volAdvMuFilter->AddNode(
            volMagnetizedFe, layer, new TGeoTranslation(0, 0, Zshift + fFeZ / 2 + layer * Fe_and_layer_z));
        if (plane == 0) {
            // Y-plane (horizontal strips along y axis)
            volAdvMuFilter->AddNode(
                ActiveLayer, layer, new TGeoTranslation(0, -layer_DY, Zshift + fFeZ + layer * Fe_and_layer_z));
        } else if (plane == 1) {
            // X-plane (vertical strips along x axis)
            volAdvMuFilter->AddNode(
                ActiveLayer,
                layer,
                new TGeoCombiTrans(TGeoTranslation(layer_DY, 0, Zshift + fFeZ + layer * Fe_and_layer_z),
                                   TGeoRotation("y_rot", 0, 0, 90)));
        }
    }   // layers
}
Bool_t AdvMuFilter::ProcessHits(FairVolume *vol)
{
    /** This method is called from the MC stepping */
    // Set parameters at entrance of volume. Reset ELoss.
    if (gMC->IsTrackEntering()) {
        fELoss = 0.;
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fPos);
        gMC->TrackMomentum(fMom);
    }
    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();

    // Create AdvMuFilterPoint at exit of active volume
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) {
        if (fELoss == 0.) {
            return kFALSE;
        }

        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();

        TParticle *p = gMC->GetStack()->GetCurrentTrack();
        Int_t pdgCode = p->GetPdgCode();
        TLorentzVector Pos;
        gMC->TrackPosition(Pos);
        TLorentzVector Mom;
        gMC->TrackMomentum(Mom);
        Int_t detID = 0;
        gMC->CurrentVolID(detID);
        const char *name;
        name = gMC->CurrentVolName();
        // Check which volume is actually hit and what detID is given
        LOG(DEBUG) << "AdvMuFilterPoint DetID " << detID << " Hit volume name " << name;
        fVolumeID = detID;
        Double_t xmean = (fPos.X() + Pos.X()) / 2.;
        Double_t ymean = (fPos.Y() + Pos.Y()) / 2.;
        Double_t zmean = (fPos.Z() + Pos.Z()) / 2.;
        AddHit(fTrackID,
               detID,
               TVector3(xmean, ymean, zmean),
               TVector3(fMom.Px(), fMom.Py(), fMom.Pz()),
               fTime,
               fLength,
               fELoss,
               pdgCode);

        // Increment number of det points in TParticle
        ShipStack *stack = (ShipStack *)gMC->GetStack();
        stack->AddPoint(kAdvSNDMuFilter);
    }
    return kTRUE;
}

void AdvMuFilter::GetPosition(Int_t detID, TVector3 &A, TVector3 &B)
{
    // Calculate the detector id as per the geofile, where strips are disrespected
    int strip = (detID) % 1024;                // actual strip ID
    int geofile_detID = detID - strip + 999;   // the det id number needed to read the geometry

    int layer = geofile_detID >> 17;
    int row = (geofile_detID >> 13) % 8;
    int column = (geofile_detID >> 11) % 4;
    int sensor = geofile_detID;
    int sensor_module = advsnd::hcal::columns * row + 1 + column;

    double local_pos[3] = {0, 0, 0};
    TString path = TString::Format("/cave_1/"
                                   "Detector_0/"
                                   "volAdvMuFilter_0/"
                                   "HCAL_Layer_%d/"
                                   "SensorModule_%d/"
                                   "HCAL_SensorVolume_%d",
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

void AdvMuFilter::EndOfEvent() { fAdvMuFilterPointCollection->Clear(); }

void AdvMuFilter::Register()
{
    /** This will create a branch in the output tree called
     AdvMuFilterPoint, setting the last parameter to kFALSE means:
     this collection will not be written to the file, it will exist
     only during the simulation.
     */

    FairRootManager::Instance()->Register("AdvMuFilterPoint", "AdvMuFilter", fAdvMuFilterPointCollection, kTRUE);
}
TClonesArray *AdvMuFilter::GetCollection(Int_t iColl) const
{
    if (iColl == 0) {
        return fAdvMuFilterPointCollection;
    } else {
        return NULL;
    }
}

void AdvMuFilter::Reset() { fAdvMuFilterPointCollection->Clear(); }

AdvMuFilterPoint *AdvMuFilter::AddHit(Int_t trackID,
                                      Int_t detID,
                                      TVector3 pos,
                                      TVector3 mom,
                                      Double_t time,
                                      Double_t length,
                                      Double_t eLoss,
                                      Int_t pdgCode)
{
    TClonesArray &clref = *fAdvMuFilterPointCollection;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) AdvMuFilterPoint(trackID, detID, pos, mom, time, length, eLoss, pdgCode);
}
