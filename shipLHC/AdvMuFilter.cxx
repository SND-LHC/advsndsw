//
//  AdvMuFilter.cxx
//
//  by D. Centanni
//  Sept 2022
//

#include "AdvMuFilter.h"
#include "AdvMuFilterPoint.h"

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

using std::cout;
using std::endl;
using std::to_string;
using std::string;

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

	/*
    Double_t fWallX = conf_floats["AdvMuFilter/WallX"];
    Double_t fWallY = conf_floats["AdvMuFilter/WallY"];
    Double_t fWallZ = conf_floats["AdvMuFilter/WallZ"];
    Double_t fPlaneX = conf_floats["AdvMuFilter/PlaneX"];
    Double_t fPlaneY = conf_floats["AdvMuFilter/PlaneY"];
    Double_t fPlaneZ = conf_floats["AdvMuFilter/PlaneZ"];
    Int_t fnPlanes = conf_ints["AdvMuFilter/nPlanes"];

    TGeoBBox *MFWall = new TGeoBBox("MFWall", fWallX/2., fWallY/2., fWallZ/2.);
    TGeoBBox *MFPlane = new TGeoBBox("MFPlane", fPlaneX/2., fPlaneY/2., fPlaneZ/2.);

    TGeoVolume *volMFWall = new TGeoVolume("volMFWall", MFWall, Fe);
    TGeoVolume *volMFPlane = new TGeoVolume("volMFPlane", MFPlane, Scint);

    volMFWall->SetLineColor(kGreen-4);
    volMFPlane->SetLineColor(kGray);
    AddSensitiveVolume(volMFPlane);

     
    TGeoVolumeAssembly *volAdvMuFilter  = new TGeoVolumeAssembly("volAdvMuFilter");

    detector->AddNode(volAdvMuFilter,1,new TGeoTranslation(-2.4244059999999976,0,0));

    // Positioning calculations, to be deleted once the whole AdvSND apparatus
    // will be positioned correctly
    TVector3 EmWall0_survey(5.35+42.2/2., 17.2+42.2/2., 288.92+10/2.+100); // cm
    Double_t TargetDiff = 100. - 63.941980;
    Double_t fTargetWallX = 100;
    // adding walls & trackers
    LOG(INFO) << " nPlanes: " << fnPlanes;

    for(int i=0; i<fnPlanes; i++)
    {
        volAdvMuFilter->AddNode(volMFWall, i, new TGeoTranslation(-EmWall0_survey.X()+(fTargetWallX-42.2)/2., EmWall0_survey.Y(), -TargetDiff+EmWall0_survey.Z()+i*(fWallZ+fPlaneZ)+fWallZ/2.));
        if(i!=(fnPlanes-1)) volAdvMuFilter->AddNode(volMFPlane, i, new TGeoTranslation(-EmWall0_survey.X()+(fTargetWallX-42.2)/2., EmWall0_survey.Y(), -TargetDiff+EmWall0_survey.Z()+i*(fWallZ+fPlaneZ)+fWallZ+fPlaneZ/2.));
    }
    LOG(INFO) <<"  MuFilter X: "<< -EmWall0_survey.X()+(fTargetWallX-42.2)/2.<<"  Y: "<< EmWall0_survey.Y()<< "  Z: "<< -TargetDiff+EmWall0_survey.Z()+fWallZ/2.;
}
	*/
	
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
  Int_t    fNBars         = conf_ints["AdvMuFilter/NBars"];
  Double_t fBarGap        = conf_floats["AdvMuFilter/BarGap"]; 
  Double_t fBarY          = (fMuonSysPlaneY-(fNBars-1)*fBarGap)/fNBars;
  
  TGeoVolumeAssembly *volAdvMuFilter  = new TGeoVolumeAssembly("volAdvMuFilter");
  
  // Moving the detector forward in the Z direction adding:
  //	-	354.362 which is the last target detector plane midpoint position
  // 	-	3 is the target detector plane half Z-dimension
  //	-	1 is the clearance between volAdvTarget and volAdvMuilter
  
  Double_t fTargetWallX = 50.; //cm
  TVector3 EmWall0_survey(5.35+42.2/2.-(fTargetWallX-42.2)/2., 17.2+42.2/2., 288.92+10/2.+100); // cm

  detector->AddNode(volAdvMuFilter,0,new TGeoTranslation(-2.4244059999999976-EmWall0_survey.X(), 38.3, 354.862+3+1+fFeZ/2.-41.895793+1.-3.854227000000008)); // hardcoded, try to find and elegant solution
  
  TGeoBBox *FeWall = new TGeoBBox("FeWall", fFeX/2., fFeY/2., fFeZ/2.);
  TGeoBBox *MuonSysPlane = new TGeoBBox("MuonSysPlane", fMuonSysPlaneX/2., fMuonSysPlaneY/2., fFeZ/2.+0.001);
  TGeoBBox *CoilSpace = new TGeoBBox("CoilSpace", fCoilX/2., fCoilY/2.+0.005, fFeZ/2.+0.05);
  TGeoBBox *Coil = new TGeoBBox("Coil", fCoilX/2., fCoilY/2., fCoilZ/2.);
  //TGeoBBox *MuonSysDet = new TGeoBBox("MuonSysDet", fMuonSysPlaneX/2., fMuonSysPlaneY/2., fFeGap/2.); //unused

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

  TGeoCompositeShape *MuonSysFe = new TGeoCompositeShape("MuonSysFe", "FeWall-MuonSysPlane-(CoilSpace:CoilUp)-(CoilSpace:CoilDown)-(FeCut:CutUpRight)-(FeCut:CutDownRight)-(FeCut:CutDownLeft)-(FeCut:CutUpLeft)");
  TGeoVolume *volFeWall = new TGeoVolume("volFeWall", MuonSysFe, Fe);
  TGeoVolume *volMagFe = new TGeoVolume("volMagFe", MuonSysPlane, Fe);
  volFeWall->SetLineColor(kGreen-4);
  volMagFe->SetLineColor(kGreen);
  
  //Double_t fField = 1.5; // Tesla
  //fField = fField/10; // kGauss (complying with GEANT3)
  Double_t fField = conf_floats["AdvMuFilter/Field"];
  LOG(INFO) << " Mag field: " << fField/10. << " Tesla" << endl;
  TGeoUniformMagField *magField = new TGeoUniformMagField(-fField,0, 0);
  TGeoGlobalMagField::Instance()->SetField(magField);
  volMagFe->SetField(magField);

  TGeoVolume *volCoil = new TGeoVolume("volCoil", Coil, Cu);
  volCoil->SetLineColor(kOrange+1);
  TGeoVolume *volVertCoil = new TGeoVolume("volVertCoil", VertCoil, Cu);
  volVertCoil->SetLineColor(kOrange+1);

  TGeoVolume *volBar = gGeoManager->MakeBox("volBar", Scint, fMuonSysPlaneX/2., fBarY/2., fFeGap/2.);
  volBar->SetLineColor(kGray-2);
  AddSensitiveVolume(volBar);
  TGeoVolume *volMuonSysDet;

  volAdvMuFilter->AddNode(volVertCoil, 0, new TGeoTranslation(0, 0, 0));
    for(int i = 0; i<fNplanes; i++)
    {
      volAdvMuFilter->AddNode(volFeWall, i, new TGeoTranslation(0, 0, (fCoilY+fFeZ)/2+i*(fFeZ+fFeGap)));
      volAdvMuFilter->AddNode(volMagFe, i, new TGeoTranslation(0, 0, (fCoilY+fFeZ)/2+i*(fFeZ+fFeGap)));
      if (i == fNplanes-1) continue;
      string name = "volMuonSysDet_"+to_string(i);
      volMuonSysDet = new TGeoVolumeAssembly(name.c_str());
      for (Int_t ibar=0; ibar< fNBars; ibar++){
        Double_t dy_bar = (fBarY+fBarGap)*ibar-fMuonSysPlaneY/2.+fBarY/2;
        volMuonSysDet->AddNode(volBar, 100*i+ibar, new TGeoTranslation(0, dy_bar, 0));
      }
      Double_t XY_angle = 0.;
      if (i%2 != 0) XY_angle = 90.;
      volAdvMuFilter->AddNode(volMuonSysDet, i, new TGeoCombiTrans( TGeoTranslation(0, 0, (fCoilY+fFeZ)/2+(fFeZ+fFeGap)/2.+i*(fFeZ+fFeGap)), TGeoRotation("Bar_rot", 0, 0, XY_angle)));
    }
  volAdvMuFilter->AddNode(volCoil, 0, new TGeoTranslation(0, fMuonSysPlaneY/2.+fCoilY/2., (fCoilY+fFeZ)/2+fCoilZ/2.-fFeZ/2.));
  volAdvMuFilter->AddNode(volCoil, 1, new TGeoTranslation(0, -fMuonSysPlaneY/2.-fCoilY/2., (fCoilY+fFeZ)/2+fCoilZ/2.-fFeZ/2.));

  volAdvMuFilter->AddNode(volVertCoil, 1, new TGeoTranslation(0, 0, (fNplanes-1)*(fFeZ+fFeGap)+fCoilY+fFeZ));
  
  // Now adding the second part of the spectrometer, TO DO: think about how to implement it in the Magnet.cxx class with all of the shapes.
  // for now just add a 1000 offset to identify Magnet hits from MuFilter ones
  Int_t    fDetIDOffset = 10000;
  //Double_t fNplanes2 = 20; // unused
  Double_t fMagnetsGap = 90; // cm
  Double_t FirstMagZ = (fNplanes-1)*(fFeZ+fFeGap)+fCoilY+fFeZ+fCoilY/2+0.01;
  //Double_t fShortCoilZ = fNplanes2*fFeZ; // unused
  Double_t fSpacing = 8.75; // cm

  Double_t fDownCutOffset     = conf_floats["AdvMuFilter/DownCutOffset"];
  Double_t fDownFeX           = conf_floats["AdvMuFilter/DownFeX"];
  Double_t fDownFeY           = conf_floats["AdvMuFilter/DownFeY"];
  Double_t fDownFeZ           = conf_floats["AdvMuFilter/DownFeZ"];
  Double_t CoilThickY         = conf_floats["AdvMuFilter/CoilThickY"];
  Double_t fDownFeYokeX       = conf_floats["AdvMuFilter/DownFeYokeX"];
  Double_t fDownFeYokeY       = conf_floats["AdvMuFilter/DownFeYokeY"];
  Double_t fDownFeCutX        = conf_floats["AdvMuFilter/DownFeCutX"];
  Double_t fDownFeCutY        = conf_floats["AdvMuFilter/DownFeCutY"];

  Double_t fIronCoreZ         = conf_floats["AdvMuFilter/IronCoreZ"];
  Double_t fIronCoreX1        = conf_floats["AdvMuFilter/IronCoreX1"];
  Double_t fIronCoreX2        = conf_floats["AdvMuFilter/IronCoreX2"];
  Double_t fIronCoreY1        = conf_floats["AdvMuFilter/IronCoreY1"];
  Double_t fIronCoreY2        = conf_floats["AdvMuFilter/IronCoreY2"];

  /* TGeoBBox *ShortCoil = new TGeoBBox("ShortCoil", fCoilX/2., fCoilY/2., fShortCoilZ/2.);
  TGeoVolume *volShortCoil = new TGeoVolume("volShortCoil", ShortCoil, Cu);
  volShortCoil->SetLineColor(kOrange+1);


  for(int i = 0; i< 20; i++)
  {
      volAdvMuFilter->AddNode(volFeWall, i+fNplanes, new TGeoTranslation(0, 0, FirstMagZ+fMagnetsGap+i*fFeZ-fFeGap+fSpacing));
      volAdvMuFilter->AddNode(volMagFe, i+fNplanes, new TGeoTranslation(0, 0, FirstMagZ+fMagnetsGap+i*fFeZ-fFeGap+fSpacing));
  }
  volAdvMuFilter->AddNode(volShortCoil, 0, new TGeoTranslation(0, fMuonSysPlaneY/2.+fCoilY/2., fShortCoilZ/2.-fFeZ/2.+fMagnetsGap+FirstMagZ-fFeGap+fSpacing));
  volAdvMuFilter->AddNode(volShortCoil, 1, new TGeoTranslation(0, -fMuonSysPlaneY/2.-fCoilY/2., fShortCoilZ/2.-fFeZ/2.+fMagnetsGap+FirstMagZ-fFeGap+fSpacing)); */

  TGeoBBox *_IronYoke = new TGeoBBox("_IronYoke", fDownFeX/2., fDownFeY/2., fDownFeZ/2.);

  Double_t _cutvers[8][2];
  _cutvers[0][0] = -fIronCoreX1/2.;
  _cutvers[0][1] = fIronCoreY1/2.;
  _cutvers[1][0] = fIronCoreX1/2.;
  _cutvers[1][1] = fIronCoreY1/2.;
  _cutvers[2][0] = fIronCoreX1/2.;
  _cutvers[2][1] = -fIronCoreY1/2.;
  _cutvers[3][0] = -fIronCoreX1/2.;
  _cutvers[3][1] = -fIronCoreY1/2.;
  _cutvers[4][0] = -fIronCoreX2/2.;
  _cutvers[4][1] = fIronCoreY2/2.;
  _cutvers[5][0] = fIronCoreX2/2.;
  _cutvers[5][1] = fIronCoreY2/2.;
  _cutvers[6][0] = fIronCoreX2/2.;
  _cutvers[6][1] = -fIronCoreY2/2.;
  _cutvers[7][0] = -fIronCoreX2/2.;
  _cutvers[7][1] = -fIronCoreY2/2.;
  TGeoArb8 *IronCore = new TGeoArb8("IronCore", fIronCoreZ/2., (Double_t *)_cutvers);
  
  Double_t CoilVers[8][2];
  CoilVers[0][0] = -fIronCoreX1/2.;
  CoilVers[0][1] = (fIronCoreY1+CoilThickY)/2.;
  CoilVers[1][0] = fIronCoreX1/2.;
  CoilVers[1][1] = (fIronCoreY1+CoilThickY)/2.;
  CoilVers[2][0] = fIronCoreX1/2.;
  CoilVers[2][1] = -(fIronCoreY1+CoilThickY)/2.;
  CoilVers[3][0] = -fIronCoreX1/2.;
  CoilVers[3][1] = -(fIronCoreY1+CoilThickY)/2.;
  CoilVers[4][0] = -fIronCoreX1/2.;
  CoilVers[4][1] = (fIronCoreY2+CoilThickY)/2.;
  CoilVers[5][0] = fIronCoreX1/2.;
  CoilVers[5][1] = (fIronCoreY2+CoilThickY)/2.;
  CoilVers[6][0] = fIronCoreX1/2.;
  CoilVers[6][1] = -(fIronCoreY2+CoilThickY)/2.;
  CoilVers[7][0] = -fIronCoreX1/2.;
  CoilVers[7][1] = -(fIronCoreY2+CoilThickY)/2.;
  TGeoArb8 *_DownCoil = new TGeoArb8("_DownCoil", fIronCoreZ/2.-0.001, (Double_t *)CoilVers);
  
  
  Double_t YokeHole[8][2];
  YokeHole[0][0] = -(fIronCoreX1+0.0001)/2.;
  YokeHole[0][1] = (fIronCoreY1+CoilThickY+0.0001)/2.;
  YokeHole[1][0] = (fIronCoreX1+0.0001)/2.;
  YokeHole[1][1] = (fIronCoreY1+CoilThickY+0.0001)/2.;
  YokeHole[2][0] = (fIronCoreX1+0.0001)/2.;
  YokeHole[2][1] = -(fIronCoreY1+CoilThickY+0.0001)/2.;
  YokeHole[3][0] = -(fIronCoreX1+0.0001)/2.;
  YokeHole[3][1] = -(fIronCoreY1+CoilThickY+0.0001)/2.;
  YokeHole[4][0] = -(fIronCoreX2+0.0001)/2.;
  YokeHole[4][1] = (fIronCoreY2+CoilThickY+0.0001)/2.;
  YokeHole[5][0] = (fIronCoreX2+0.0001)/2.;
  YokeHole[5][1] = (fIronCoreY2+CoilThickY+0.0001)/2.;
  YokeHole[6][0] = (fIronCoreX2+0.0001)/2.;
  YokeHole[6][1] = -(fIronCoreY2+CoilThickY+0.0001)/2.;
  YokeHole[7][0] = -(fIronCoreX2+0.0001)/2.;
  YokeHole[7][1] = -(fIronCoreY2+CoilThickY+0.0001)/2.;
  TGeoArb8 *_YokeHole = new TGeoArb8("_YokeHole", fIronCoreZ/2.+0.001, (Double_t *)YokeHole);

  Double_t downcutvers[8][2];
  downcutvers[0][0] = 0;
  downcutvers[0][1] = 0;
  downcutvers[1][0] = 0;
  downcutvers[1][1] = -fDownFeCutY;
  downcutvers[2][0] = 0;
  downcutvers[2][1] = -fDownFeCutY;
  downcutvers[3][0] = +fDownFeCutX;
  downcutvers[3][1] = 0;

  downcutvers[4][0] = 0;
  downcutvers[4][1] = 0;
  downcutvers[5][0] = 0;
  downcutvers[5][1] = -fDownFeCutY;
  downcutvers[6][0] = 0;
  downcutvers[6][1] = -fDownFeCutY;
  downcutvers[7][0] = +fDownFeCutX;
  downcutvers[7][1] = 0;
  TGeoArb8 *DownFeCut = new TGeoArb8("DownFeCut", fDownFeZ/2.+0.001, (Double_t *)downcutvers);

    TGeoTranslation *DownCutUpRight = new TGeoTranslation("DownCutUpRight", -fDownFeX/2.-0.001, fDownFeY/2.+0.001, 0);
    DownCutUpRight->RegisterYourself();
    TGeoCombiTrans *DownCutDownRight = new TGeoCombiTrans( TGeoTranslation(-fDownFeX/2.-0.001, -fDownFeY/2.-0.001, 0), TGeoRotation("rot", 0, 0, 90));
    DownCutDownRight->SetName("DownCutDownRight");
    DownCutDownRight->RegisterYourself();
    TGeoCombiTrans *DownCutDownLeft = new TGeoCombiTrans(TGeoTranslation(+fDownFeX/2.+0.001, -fDownFeY/2.-0.001, 0), TGeoRotation("rot1", 0, 0, 180));
    DownCutDownLeft->SetName("DownCutDownLeft");
    DownCutDownLeft->RegisterYourself();
    TGeoCombiTrans *DownCutUpLeft = new TGeoCombiTrans(TGeoTranslation(+fDownFeX/2.+0.001, +fDownFeY/2.+0.001, 0), TGeoRotation("rot2", 0, 0, -90));
    DownCutUpLeft->SetName("DownCutUpLeft");
    DownCutUpLeft->RegisterYourself();
  
  TGeoVolume *volIronCore = new TGeoVolume("volIronCore", IronCore, Fe);

  TGeoCompositeShape *DownCoil = new TGeoCompositeShape("DownCoil", "_DownCoil-IronCore");
  TGeoVolume *volDownCoil = new TGeoVolume("volDownCoil", DownCoil, Cu);


  TGeoCompositeShape *IronYoke = new TGeoCompositeShape("IronYoke", "_IronYoke-_YokeHole-(DownFeCut:DownCutUpRight)-(DownFeCut:DownCutDownRight)-(DownFeCut:DownCutDownLeft)-(DownFeCut:DownCutUpLeft)");

  
  TGeoVolume *volIronYoke = new TGeoVolume("volIronYoke", IronYoke, Fe);

  Double_t fDownField = conf_floats["AdvMuFilter/DownField"];
  LOG(INFO) << " Downstream Mag field: " << fDownField/10. << " Tesla" << endl;
  //TGeoUniformMagField *DownmagField = new TGeoUniformMagField(-fDownField,0, 0);
  //TGeoGlobalMagField::Instance()->SetField(DownmagField);
  volIronCore->SetField(magField);

  volIronYoke->SetLineColor(kGreen-4);
  volIronYoke->SetTransparency(60);
  volIronCore->SetLineColor(kGreen);

  TGeoBBox *DownVertCoil1 = new TGeoBBox("DownVertCoil1", fIronCoreX1/2., fIronCoreY1/2., CoilThickY/4.);
  TGeoBBox *DownVertCoil2 = new TGeoBBox("DownVertCoil2", fIronCoreX2/2., fIronCoreY2/2., CoilThickY/4.);

  TGeoVolume *volDownVertCoil1 = new TGeoVolume("volDownVertCoil1", DownVertCoil1, Cu);
  TGeoVolume *volDownVertCoil2 = new TGeoVolume("volDownVertCoil2", DownVertCoil2, Cu);

  volDownCoil->SetLineColor(kOrange+1);
  volDownCoil->SetLineColor(kOrange+1);
  volDownVertCoil1->SetLineColor(kOrange+1);
  volDownVertCoil2->SetLineColor(kOrange+1);

  TGeoVolumeAssembly *volDownstreamMagnet  = new TGeoVolumeAssembly("volDownstreamMagnet");
  volDownstreamMagnet->AddNode(volDownVertCoil1, 0, new TGeoTranslation(0, 0, -fDownFeZ/2.-CoilThickY/4));
  volDownstreamMagnet->AddNode(volDownVertCoil2, 0, new TGeoTranslation(0, 0, fDownFeZ/2.+CoilThickY/4));
  volDownstreamMagnet->AddNode(volIronYoke, 0, 0);
  volDownstreamMagnet->AddNode(volIronCore, 0, 0);
  volDownstreamMagnet->AddNode(volDownCoil, 0, 0);


  // Trackers part
  Double_t fMagTrackerZ = conf_floats["AdvMuFilter/MagTrackerZ"]; // cm Alu tubes diameter
  TGeoBBox *MagTracker = new TGeoBBox("MagTracker", 80/2., 80/2., fMagTrackerZ/2.);
  //TGeoBBox *MagTracker = new TGeoBBox("MagTracker", fMuonSysPlaneX/2., fMuonSysPlaneY/2., fMagTrackerZ/2.);
  TGeoVolume *volMagTracker = new TGeoVolume("volMagTracker", MagTracker, Al);
  volMagTracker->SetLineColor(kGray);
  AddSensitiveVolume(volMagTracker);
  TGeoBBox *DownMagTracker = new TGeoBBox("DownMagTracker", fIronCoreX2/2., fIronCoreY2/2., fMagTrackerZ/2.);
  TGeoVolume *volDownMagTracker = new TGeoVolume("volDownMagTracker", DownMagTracker, Al);
  volDownMagTracker->SetLineColor(kGray);
  AddSensitiveVolume(volDownMagTracker);

  volAdvMuFilter->AddNode(volMagTracker, 0+fDetIDOffset, new TGeoTranslation(0, 0, FirstMagZ+fMagTrackerZ/2));
  volAdvMuFilter->AddNode(volMagTracker, 1+fDetIDOffset, new TGeoTranslation(0, 0, FirstMagZ+fMagnetsGap+fMagTrackerZ/2)); // Somehow magtrackers are not rendered in OGL, hardcoding the displacement
  volAdvMuFilter->AddNode(volDownstreamMagnet, 0, new TGeoTranslation(0, 0, FirstMagZ+fMagnetsGap+fMagTrackerZ+CoilThickY/2+fDownFeZ/2.+0.01+2*0.87));
  volAdvMuFilter->AddNode(volDownMagTracker, 2+fDetIDOffset, new TGeoTranslation(0, 0, FirstMagZ+fMagnetsGap+fMagTrackerZ+CoilThickY+fDownFeZ+fMagTrackerZ/2+0.02+2*1.74));

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
