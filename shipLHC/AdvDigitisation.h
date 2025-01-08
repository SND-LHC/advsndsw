#ifndef SHIPLHC_ADVDIGITISATION_H_
#define SHIPLHC_ADVDIGITISATION_H_

#include "AdvTargetPoint.h"
#include "AdvSignal.h"
#include "TVector3.h"
#include "TGeoNavigator.h"

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

#include <iostream>
#include <vector>
using namespace std;



class AdvDigitisation
{
  public:
    AdvDigitisation();
    void digirun(Int_t detID, const std::vector<AdvTargetPoint*>& V);
    TVector3 getLocal(Int_t detID, TVector3 global_pos);
    std::vector<AdvSignal> digirunoutput(Int_t detID, const std::vector<AdvTargetPoint*>& V);
    
};
#endif
