//
//  AdvMuFilter.cxx
//
//  by D. Centanni and O. Lantwin
//  2024
//

#include "AdvMuFilter.h"
#include "AdvMuFilterPoint.h"
#include "SiSensor.h"

#include "TGeoManager.h"
#include "FairRun.h"                    // for FairRun
#include "FairRuntimeDb.h"              // for FairRuntimeDb
#include <iosfwd>                    // for ostream
#include "TList.h"                      // for TListIter, TList (ptr only)
#include "TObjArray.h"                  // for TObjArray
#include "TString.h"                    // for TString

#include "TGeoBBox.h"
#include "TGeoTrd1.h"
#include "TGeoSphere.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoTrd1.h"
#include "TGeoArb8.h"
#include "TParticle.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoTransform.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include "ShipDetectorList.h"
#include "ShipUnit.h"
#include "ShipStack.h"

#include "TGeoTrd2.h" 
#include "TGeoCompositeShape.h"

#include "TGeoUniformMagField.h"
#include "TGeoGlobalMagField.h"
#include "TVector3.h"
#include <stddef.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream,etc
#include <string.h>
#include <cstring>
#include <ROOT/TSeq.hxx>

using std::cout;
using std::endl;
using std::to_string;
using std::string;
using ROOT::TSeq;

using namespace ShipUnit;

AdvMuFilter::AdvMuFilter()
: FairDetector("AdvMuFilter", "",kTRUE),
  fTrackID(-1),
  fVolumeID(-1),
  fPos(),
  fMom(),
  fTime(-1.),
  fLength(-1.),
  fELoss(-1),
  fAdvMuFilterPointCollection(new TClonesArray("AdvMuFilterPoint"))
{
}

AdvMuFilter::AdvMuFilter(const char* name, Bool_t Active,const char* Title)
: FairDetector(name, true, kAdvSNDMuFilter),
  fTrackID(-1),
  fVolumeID(-1),
  fPos(),
  fMom(),
  fTime(-1.),
  fLength(-1.),
  fELoss(-1),
  fAdvMuFilterPointCollection(new TClonesArray("AdvMuFilterPoint"))
{
}

AdvMuFilter::~AdvMuFilter()
{
    if (fAdvMuFilterPointCollection) {
        fAdvMuFilterPointCollection->Delete();
        delete fAdvMuFilterPointCollection;
    }
}

void AdvMuFilter::Initialize()
{
  FairDetector::Initialize();
}


// -----   Private method InitMedium
Int_t AdvMuFilter::InitMedium(const char* name)
{
  static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
  static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
  static FairGeoMedia *media=geoFace->getMedia();
  static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();
    
  FairGeoMedium *ShipMedium=media->getMedium(name);
    
  if (!ShipMedium)
    {
      Fatal("InitMedium","Material %s not defined in media file.", name);
      return -1111;
    }
  TGeoMedium* medium=gGeoManager->GetMedium(name);
  if (medium!=NULL)
    return ShipMedium->getMediumIndex();
  return geoBuild->createMedium(ShipMedium);
}

void AdvMuFilter::ConstructGeometry()
{
    // Geometry implementation from D. Centanni
    TGeoVolume *top=gGeoManager->GetTopVolume();
    TGeoVolume *detector = gGeoManager->FindVolumeFast("Detector");
    if(!detector)  LOG(ERROR) << "no Detector volume found " ;

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
	TGeoMedium *Al =gGeoManager->GetMedium("aluminium");

  InitMedium("silicon");
  TGeoMedium *Silicon = gGeoManager->GetMedium("silicon");

	Double_t fMuonSysPlaneX = conf_floats["AdvMuFilter/MuonSysPlaneX"];
  Double_t fMuonSysPlaneY = conf_floats["AdvMuFilter/MuonSysPlaneY"]; 
  Double_t fCutOffset     = conf_floats["AdvMuFilter/CutOffset"];  
  Double_t fFeX           = conf_floats["AdvMuFilter/FeX"];
  Double_t fFeY           = conf_floats["AdvMuFilter/FeY"];
  Double_t fFeZ           = conf_floats["AdvMuFilter/FeZ"];  
  Double_t fFeGap         = conf_floats["AdvMuFilter/FeGap"];  
  Int_t    fNplanes       = conf_ints["AdvMuFilter/Nplanes"];  
  Double_t fCoilX         = conf_floats["AdvMuFilter/CoilX"];
  Double_t fCoilY         = conf_floats["AdvMuFilter/CoilY"];  
  Double_t fCoilZ         = (fNplanes)*(fFeZ+fFeGap)-fFeGap;
  Double_t fFeYokeX       = (fFeX-fMuonSysPlaneX)/2.;
  Double_t fFeYokeY       = (fFeY-fMuonSysPlaneY-fCoilY)/2.;
  Double_t fFeCutX        = fFeYokeX - fCutOffset;
  Double_t fFeCutY        = fFeYokeY - fCutOffset;

  
  TGeoVolumeAssembly *volAdvMuFilter  = new TGeoVolumeAssembly("volAdvMuFilter");
  
  // Moving the detector forward in the Z direction adding:
  //	-	354.362 which is the last target detector plane midpoint position
  // 	-	3 is the target detector plane half Z-dimension
  //	-	1 is the clearance between volAdvTarget and volAdvMuilter
  
  Double_t fTargetWallX = 50.; //cm
  TVector3 EmWall0_survey(5.35+42.2/2.-(fTargetWallX-42.2)/2., 17.2+42.2/2., 288.92+10/2.+100); // cm

  detector->AddNode(volAdvMuFilter,0,new TGeoTranslation(-2.4244059999999976-EmWall0_survey.X(), 38.3, 354.862+3+1+fFeZ/2.-41.895793+1.-3.854227000000008+3.7497190000000046)); // hardcoded, try to find and elegant solution
  
  TGeoBBox *FeWall = new TGeoBBox("FeWall", fFeX/2., fFeY/2., fFeZ/2.);
  TGeoBBox *MuonSysPlane = new TGeoBBox("MuonSysPlane", fMuonSysPlaneX/2., fMuonSysPlaneY/2., fFeZ/2.+0.001);
  TGeoBBox *CoilSpace = new TGeoBBox("CoilSpace", fCoilX/2., fCoilY/2.+0.005, fFeZ/2.+0.05);
  TGeoBBox *Coil = new TGeoBBox("Coil", fCoilX/2., fCoilY/2., fCoilZ/2.);

  TGeoBBox *VertCoil = new TGeoBBox("VertCoil", fCoilX/2., fMuonSysPlaneY/2., fCoilY/2.);

  Double_t cutvers[8][2];
  cutvers[0][0] = 0;
  cutvers[0][1] = 0;
  cutvers[1][0] = 0;
  cutvers[1][1] = -fFeCutY;
  cutvers[2][0] = 0;
  cutvers[2][1] = -fFeCutY;
  cutvers[3][0] = +fFeCutX;
  cutvers[3][1] = 0;

  cutvers[4][0] = 0;
  cutvers[4][1] = 0;
  cutvers[5][0] = 0;
  cutvers[5][1] = -fFeCutY;
  cutvers[6][0] = 0;
  cutvers[6][1] = -fFeCutY;
  cutvers[7][0] = +fFeCutX;
  cutvers[7][1] = 0;
  TGeoArb8 *FeCut = new TGeoArb8("FeCut", fFeZ/2.+0.001, (Double_t *)cutvers);

  TGeoTranslation *CutUpRight = new TGeoTranslation("CutUpRight", -fFeX/2.-0.001, fFeY/2.+0.001, 0);
  CutUpRight->RegisterYourself();
  TGeoCombiTrans *CutDownRight = new TGeoCombiTrans( TGeoTranslation(-fFeX/2.-0.001, -fFeY/2.-0.001, 0), TGeoRotation("rot", 0, 0, 90));
  CutDownRight->SetName("CutDownRight");
  CutDownRight->RegisterYourself();
  TGeoCombiTrans *CutDownLeft = new TGeoCombiTrans(TGeoTranslation(+fFeX/2.+0.001, -fFeY/2.-0.001, 0), TGeoRotation("rot1", 0, 0, 180));
  CutDownLeft->SetName("CutDownLeft");
  CutDownLeft->RegisterYourself();
  TGeoCombiTrans *CutUpLeft = new TGeoCombiTrans(TGeoTranslation(+fFeX/2.+0.001, +fFeY/2.+0.001, 0), TGeoRotation("rot2", 0, 0, -90));
  CutUpLeft->SetName("CutUpLeft");
  CutUpLeft->RegisterYourself();


  TGeoTranslation *CoilUp = new TGeoTranslation("CoilUp", 0, fMuonSysPlaneY/2.+fCoilY/2., 0);
  TGeoTranslation *CoilDown = new TGeoTranslation("CoilDown", 0, -fMuonSysPlaneY/2.-fCoilY/2., 0);
  CoilUp->RegisterYourself();
  CoilDown->RegisterYourself();

  TGeoCompositeShape *MuonSysFe = new TGeoCompositeShape("MuonSysFe", "FeWall-MuonSysPlane-(CoilSpace:CoilUp)-(CoilSpace:CoilDown)");
  TGeoVolume *volFeWall = new TGeoVolume("volFeWall", MuonSysFe, Fe);
  TGeoVolume *volMagFe = new TGeoVolume("volMagFe", MuonSysPlane, Fe);
  volFeWall->SetLineColor(kGreen-4);
  volMagFe->SetLineColor(kGreen);
  
  Double_t fField = conf_floats["AdvMuFilter/Field"];
  LOG(INFO) << " Mag field: " << fField/10. << " Tesla" << endl;
  TGeoUniformMagField *magField = new TGeoUniformMagField(-fField,0, 0);
  TGeoGlobalMagField::Instance()->SetField(magField);
  volMagFe->SetField(magField);

  TGeoVolume *volCoil = new TGeoVolume("volCoil", Coil, Cu);
  volCoil->SetLineColor(kOrange+1);
  TGeoVolume *volVertCoil = new TGeoVolume("volVertCoil", VertCoil, Cu);
  volVertCoil->SetLineColor(kOrange+1);

  // Minimal configuration includes Silicon strip detectors, for now naive copy of AdvTarget layout
  // Silicon tracker module
  //
  // See https://indico.cern.ch/event/1201858/#81-detector-simulation for technical diagrams and renders
  //
  // Passive part
  TGeoBBox *Support = new TGeoBBox("Support", advsnd::module_length / 2, advsnd::module_width / 2, 3.0 * mm / 2);
  TGeoVolume *SupportVolume = new TGeoVolume("SupportVolume", Support, Polystyrene);
  SupportVolume->SetLineColor(kGray);
  // Active part
  TGeoBBox *SensorShape = new TGeoBBox("SensorShape", advsnd::sensor_length / 2, advsnd::sensor_width / 2, 0.5 * mm / 2);
  TGeoVolume *SensorVolume = new TGeoVolume("SensorVolume", SensorShape, Silicon);
  SensorVolume->SetLineColor(kGreen);
  AddSensitiveVolume(SensorVolume);

  volAdvMuFilter->AddNode(volVertCoil, 0, new TGeoTranslation(0, 0, 0));
    for(int i = 0; i<fNplanes; i++)
    {
      volAdvMuFilter->AddNode(volFeWall, i, new TGeoTranslation(0, 0, (fCoilY+fFeZ)/2+i*(fFeZ+fFeGap)));
      volAdvMuFilter->AddNode(volMagFe, i, new TGeoTranslation(0, 0, (fCoilY+fFeZ)/2+i*(fFeZ+fFeGap)));
      if (i == fNplanes-1) continue;

      int station = i;
      TGeoVolumeAssembly *TrackingStation = new TGeoVolumeAssembly("TrackingStation");

      if (i < advsnd::hcal::stations) {
        using namespace advsnd::hcal;

        // Each tracking station consists of X and Y planes
        for (auto &&plane : TSeq(planes)) {
          TGeoVolumeAssembly *TrackerPlane = new TGeoVolumeAssembly("TrackerPlane");
          int j = 0;
            // Each plane consists of 4 modules
            for (auto &&row : TSeq(rows)) {
                for (auto &&column : TSeq(columns)) {
                    // Each module in turn consists of two sensors on a support
                    TGeoVolumeAssembly *SensorModule = new TGeoVolumeAssembly("SensorModule");
                    SensorModule->AddNode(SupportVolume, 1);
                    for (auto &&sensor : TSeq(advsnd::sensors)) {
                        int sensor_id =  (station << 5) + (plane << 4) + (row << 2) + (column << 1) + sensor;
                        SensorModule->AddNode(SensorVolume,
                                              sensor_id,
                                              new TGeoTranslation(-advsnd::module_length / 2 + 46.95 * mm + advsnd::sensor_length / 2
                                                                      + sensor * (advsnd::sensor_length + advsnd::sensor_gap),
                                                                  0,
                                                                  +3 * mm / 2 + 0.5 * mm / 2));
                    }
                    TrackerPlane->AddNode(
                        SensorModule,
                        ++j,
                        new TGeoCombiTrans(
                                  // Offset modules as needed by row and column
                                  TGeoTranslation(-plane_width / 2. + advsnd::module_length / 2. + column * (advsnd::module_length + module_column_gap) - (column == 2 ? 43.85 * mm : 3.35 * mm),
                                          (row-(rows-1.)/2.) * (advsnd::module_width + module_row_gap),
                                          (column-(columns-1.)/2.) * 4 * mm),
                            // Rotate every module of the second column
                            TGeoRotation(TString::Format("rot%d", j), 0, 0, (column != 2) ? 180 : 0)));
                }
            }
            if (plane == 0) {
                // X-plane
                TrackingStation->AddNode(TrackerPlane, plane);
            } else if (plane == 1) {
                // Y-plane
                TrackingStation->AddNode(
                    TrackerPlane,
                    plane,
                    new TGeoCombiTrans(TGeoTranslation(0, 0, +3.5 * mm + 0.5 * mm), TGeoRotation("y_rot", 0, 0, 90)));
            }
        }

      }
      else if (i >= advsnd::hcal::stations) {
        using namespace advsnd::muon;

        // Each tracking station consists of X and Y planes
        for (auto &&plane : TSeq(planes)) {
            TGeoVolumeAssembly *TrackerPlane = new TGeoVolumeAssembly("TrackerPlane");
            int j = 0;
            // Each plane consists of 4 modules
            for (auto &&row : TSeq(rows)) {
                for (auto &&column : TSeq(columns)) {
                    // Each module in turn consists of two sensors on a support
                    TGeoVolumeAssembly *SensorModule = new TGeoVolumeAssembly("SensorModule");
                    SensorModule->AddNode(SupportVolume, 1);
                    for (auto &&sensor : TSeq(advsnd::sensors)) {
                        int sensor_id =  (station << 5) + (plane << 4) + (row << 2) + (column << 1) + sensor;
                        SensorModule->AddNode(SensorVolume,
                                              sensor_id,
                                              new TGeoTranslation(-advsnd::module_length / 2 + 46.95 * mm + advsnd::sensor_length / 2
                                                                      + sensor * (advsnd::sensor_length + advsnd::sensor_gap),
                                                                  0,
                                                                  +3 * mm / 2 + 0.5 * mm / 2));
                    }
                    TrackerPlane->AddNode(
                        SensorModule,
                        ++j,
                        new TGeoCombiTrans(
                            // Offset modules as needed by row and column
                          TGeoTranslation(-plane_width / 2. +advsnd::module_length / 2. + column * (advsnd::module_length + module_column_gap) - (column == 2 ? 43.85 * mm : 3.35 * mm),
                                          (row-(rows-1.)/2.) * (advsnd::module_width + module_row_gap),
                                          (column-(columns-1.)/2.) * 4 * mm),
                            // Rotate every module of the second column
                          TGeoRotation(TString::Format("rot%d", j), 0, 0, (column != 2) ? 180 : 0)));
                }
            }
            if (plane == 0) {
                // X-plane
                TrackingStation->AddNode(TrackerPlane, plane);
            } else if (plane == 1) {
                // Y-plane
                TrackingStation->AddNode(
                    TrackerPlane,
                    plane,
                    new TGeoCombiTrans(TGeoTranslation(0, 0, +3.5 * mm + 0.5 * mm), TGeoRotation("y_rot", 0, 0, 90)));
            }
        }

      }
      volAdvMuFilter->AddNode(TrackingStation, i, new TGeoTranslation(0, 0, (fCoilY+fFeZ)/2+(fFeZ+fFeGap)/2.+i*(fFeZ+fFeGap)));
    }
  volAdvMuFilter->AddNode(volCoil, 0, new TGeoTranslation(0, fMuonSysPlaneY/2.+fCoilY/2., (fCoilY+fFeZ)/2+fCoilZ/2.-fFeZ/2.));
  volAdvMuFilter->AddNode(volCoil, 1, new TGeoTranslation(0, -fMuonSysPlaneY/2.-fCoilY/2., (fCoilY+fFeZ)/2+fCoilZ/2.-fFeZ/2.));
  volAdvMuFilter->AddNode(volVertCoil, 1, new TGeoTranslation(0, 0, (fNplanes-1)*(fFeZ+fFeGap)+fCoilY+fFeZ));

}
Bool_t  AdvMuFilter::ProcessHits(FairVolume* vol)
{
  /** This method is called from the MC stepping */
  //Set parameters at entrance of volume. Reset ELoss.
  if ( gMC->IsTrackEntering() ) {
    fELoss  = 0.;
    fTime   = gMC->TrackTime() * 1.0e09;
    fLength = gMC->TrackLength();
    gMC->TrackPosition(fPos);
    gMC->TrackMomentum(fMom);
  }
  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();

  // Create AdvMuFilterPoint at exit of active volume
  if (  gMC->IsTrackExiting()    ||
        gMC->IsTrackStop()       ||
        gMC->IsTrackDisappeared()   ) {
        if (fELoss == 0. ) {return kFALSE; }

    fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();

    TParticle* p=gMC->GetStack()->GetCurrentTrack();
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
    LOG(DEBUG)<< "AdvMuFilterPoint DetID " << detID << " Hit volume name " << name;
    fVolumeID = detID;
    Double_t xmean = (fPos.X()+Pos.X())/2. ;
	  Double_t ymean = (fPos.Y()+Pos.Y())/2. ;
	  Double_t zmean = (fPos.Z()+Pos.Z())/2. ;
    AddHit(fTrackID, detID,  TVector3(xmean, ymean,  zmean),
		TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,
			fELoss, pdgCode);

    // Increment number of det points in TParticle
		ShipStack* stack = (ShipStack*) gMC->GetStack();
		stack->AddPoint(kAdvSNDMuFilter);
    }
        return kTRUE;
}

void AdvMuFilter::EndOfEvent()
{
    fAdvMuFilterPointCollection->Clear();
}

void AdvMuFilter::Register()
{
    /** This will create a branch in the output tree called
     AdvMuFilterPoint, setting the last parameter to kFALSE means:
     this collection will not be written to the file, it will exist
     only during the simulation.
     */

    FairRootManager::Instance()->Register("AdvMuFilterPoint", "AdvMuFilter",
                                            fAdvMuFilterPointCollection, kTRUE);
}
TClonesArray* AdvMuFilter::GetCollection(Int_t iColl) const
{
    if (iColl == 0) { return fAdvMuFilterPointCollection; }
    else { return NULL; }
}

void AdvMuFilter::Reset()
{
    fAdvMuFilterPointCollection->Clear();
}

AdvMuFilterPoint* AdvMuFilter::AddHit(Int_t trackID, Int_t detID,
                           TVector3 pos, TVector3 mom,
                           Double_t time, Double_t length,
                           Double_t eLoss, Int_t pdgCode)
{
    TClonesArray& clref = *fAdvMuFilterPointCollection;
    Int_t size = clref.GetEntriesFast();
    return new(clref[size]) AdvMuFilterPoint(trackID, detID, pos, mom,
                                       time, length, eLoss, pdgCode);
}


ClassImp(AdvMuFilter)
