//
//  AdvMuFilter.cxx
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

    InitMedium("polyvinyltoluene");
    TGeoMedium *Scint = gGeoManager->GetMedium("polyvinyltoluene");

    InitMedium("CoilCopper");
    TGeoMedium *Cu = gGeoManager->GetMedium("CoilCopper");

    InitMedium("Polystyrene");
    TGeoMedium *Polystyrene = gGeoManager->GetMedium("Polystyrene");

    InitMedium("aluminium");
    TGeoMedium *Al = gGeoManager->GetMedium("aluminium");

    InitMedium("silicon");
    TGeoMedium *Silicon = gGeoManager->GetMedium("silicon");

    Double_t fMuonSysPlaneX = conf_floats["AdvMuFilter/MuonSysPlaneX"];
    Double_t fMuonSysPlaneY = conf_floats["AdvMuFilter/MuonSysPlaneY"];
    Double_t fMuonSysPlaneZ = conf_floats["AdvMuFilter/MuonSysPlaneZ"];
    Double_t fFeX = conf_floats["AdvMuFilter/FeX"];
    Double_t fFeY = conf_floats["AdvMuFilter/FeY"];
    Double_t fFeZ = conf_floats["AdvMuFilter/FeZ"];
    Double_t fFeGap = conf_floats["AdvMuFilter/FeGap"];
    Double_t fCurvRadius = conf_floats["AdvMuFilter/CurvRadius"];
    Int_t fNplanes = conf_ints["AdvMuFilter/Nplanes"];
    Double_t fCoilZ = conf_floats["AdvMuFilter/CoilZ"];
    Double_t fCoilY = conf_floats["AdvMuFilter/CoilY"];
    Double_t fCoilX = fMuonSysPlaneX+2*fFeGap-2*(fCoilZ+fCurvRadius);
    Double_t fFeYokeX = (fFeX-fMuonSysPlaneX-2*fFeGap)/2.;
    Double_t fFeYokeY = (fFeY-fMuonSysPlaneY)/2.;


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
                      new TGeoTranslation(-2.4244059999999976 - EmWall0_survey.X(),
                                          38.3,
                                          354.862 + 3 + 1 + fFeZ / 2. - 41.895793 + 1. - 3.854227000000008
                                              + 3.7497190000000046 - 64.414875));   // hardcoded, try to find and elegant solution
    
    // Iron and Detector's shapes
    TGeoBBox *FeBlock = new TGeoBBox("FeBlock", fFeX/2., fFeY/2., fFeZ/2.);
    TGeoBBox *FeGap   = new TGeoBBox("FeGap", (2*fFeGap+fMuonSysPlaneX)/2, fMuonSysPlaneY/2., fFeZ/2.+0.001);
    TGeoCompositeShape *FeSlab = new TGeoCompositeShape("FeSlab", "FeBlock-FeGap");
    TGeoBBox *MagFe   = new TGeoBBox("MagFe", fMuonSysPlaneX/2., fMuonSysPlaneY/2., fFeZ/2.);
    TGeoBBox *MuonSysPlane   = new TGeoBBox("MuonSysPlane", fMuonSysPlaneX/2., fMuonSysPlaneY/2., fMuonSysPlaneZ/2.);

    // Coil shapes
    TGeoBBox *FrontCoil = new TGeoBBox("FrontCoil", fCoilX/2., fCoilY/2., fCoilZ/2.);
    TGeoBBox *LatCoil = new TGeoBBox("LatCoil", fCoilZ/2., fCoilY/2., fNplanes*(fFeZ+2*fMuonSysPlaneZ)/2.);
    TGeoTubeSeg *CurvCoil = new TGeoTubeSeg("CurvCoil", fCurvRadius, fCurvRadius+fCoilZ, fCoilY/2., 0, 90);

    TGeoCombiTrans *Left = new TGeoCombiTrans(TGeoTranslation(-fCoilX/2., 0, -(fCoilZ+fCurvRadius)/2.-fCurvRadius/2.), TGeoRotation("rot1", 0, 90, 90));
    Left->SetName("Left");
    Left->RegisterYourself();
    TGeoCombiTrans *Right = new TGeoCombiTrans(TGeoTranslation(+fCoilX/2., 0, -(fCoilZ+fCurvRadius)/2.-fCurvRadius/2.), TGeoRotation("rot2", 0, 90, 0));
    Right->SetName("Right");
    Right->RegisterYourself();
    TGeoCompositeShape *FrontCoilShape = new TGeoCompositeShape("FrontCoilShape", "FrontCoil+(CurvCoil:Left)+(CurvCoil:Right)");

    TGeoVolume *volFrontCoil = new TGeoVolume("volFrontCoil", FrontCoilShape, Cu);
    volFrontCoil->SetLineColor(kOrange+1);
    volAdvMuFilter->AddNode(volFrontCoil, 0, new TGeoCombiTrans(TGeoTranslation(0, 0, fCoilZ/2.), TGeoRotation("rot3", 0, 180, 0)));
    volAdvMuFilter->AddNode(volFrontCoil, 1, new TGeoTranslation(0, 0, fCoilZ+2*fCurvRadius+fCoilZ/2.+fNplanes*(fFeZ+2*fMuonSysPlaneZ)));
    TGeoTranslation *LatLeft = new TGeoTranslation(-fCoilX/2.-fCoilZ/2.-fCurvRadius, 0, 0);
    LatLeft->SetName("LatLeft");
    LatLeft->RegisterYourself();
    TGeoTranslation *LatRight = new TGeoTranslation(+fCoilX/2.+fCoilZ/2.+fCurvRadius, 0, 0);
    LatRight->SetName("LatRight");
    LatRight->RegisterYourself();
    TGeoCompositeShape *LatCoilShape = new TGeoCompositeShape("LatCoilShape", "(LatCoil:LatLeft)+(LatCoil:LatRight)");
    TGeoVolume *volLatCoil = new TGeoVolume("volLatCoil", LatCoilShape, Cu);
    volLatCoil->SetLineColor(kOrange+1);
    volAdvMuFilter->AddNode(volLatCoil, 0, new TGeoTranslation(0, 0, fCoilZ+fCurvRadius+(fNplanes/2)*(fFeZ+2*fMuonSysPlaneZ)));

    TGeoVolume *volMuonSysPlane = new TGeoVolume("volMuonSysPlane", MuonSysPlane, Silicon);
    AddSensitiveVolume(volMuonSysPlane);
    TGeoVolume *volFeSlab = new TGeoVolume("volFeSlab", FeSlab, Fe);
    TGeoVolume *volMagFe = new TGeoVolume("volMagFe", MagFe, Fe);

    Double_t fField = conf_floats["AdvMuFilter/Field"];
    LOG(INFO) << " Mag field: " << fField / 10. << " Tesla" << endl;
    TGeoUniformMagField *magField = new TGeoUniformMagField(0, fField, 0);
    TGeoGlobalMagField::Instance()->SetField(magField);
    volMagFe->SetField(magField);

    volMuonSysPlane->SetLineColor(kGray-2);
    volFeSlab->SetLineColor(kGreen-4);
    volMagFe->SetLineColor(kGreen);

    volAdvMuFilter->AddNode(volMuonSysPlane, 0, new TGeoTranslation(0, 0, fMuonSysPlaneZ/2+(fCoilZ+fCurvRadius)));
    for (auto plane = 1; plane<fNplanes+1; plane++) {
        Double_t station_length = 2*fMuonSysPlaneZ+fFeZ;
        Double_t Zshift = fMuonSysPlaneZ+fCoilZ+fCurvRadius;
        volAdvMuFilter->AddNode(volFeSlab, plane, new TGeoTranslation(0, 0, Zshift+fFeZ/2+(plane-1)*station_length));
        volAdvMuFilter->AddNode(volMagFe, plane, new TGeoTranslation(0, 0, Zshift+fFeZ/2+(plane-1)*station_length));
        volAdvMuFilter->AddNode(volMuonSysPlane, plane*10, new TGeoTranslation(0, 0, Zshift+fFeZ+fMuonSysPlaneZ/2+(plane-1)*station_length));
        if (plane == fNplanes) continue;
        volAdvMuFilter->AddNode(volMuonSysPlane, plane*10+1, new TGeoTranslation(0, 0, Zshift+fFeZ+fMuonSysPlaneZ+fMuonSysPlaneZ/2+(plane-1)*station_length));
    }
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

    int station = geofile_detID >> 17;
    int plane = (geofile_detID >> 16) % 2;
    int row = (geofile_detID >> 13) % 8;
    int column = (geofile_detID >> 11) % 4;
    int sensor = geofile_detID;
    int sensor_module = advsnd::muon::columns * row + 1 + column;

    double global_pos[3];
    double local_pos[3] = {0, 0, 0};
    TString path = TString::Format("/cave_1/"
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
