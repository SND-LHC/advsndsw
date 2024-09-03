#include "ChargeDivision.h"

#include "AdvTargetPoint.h"
#include "EnergyFluctUnit.h"
#include "SiG4UniversalFluctuation.h"
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


#include <TDatabasePDG.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

ChargeDivision::ChargeDivision() {}

void ChargeDivision::ReadPulseShape(std::string PulseFileName)
{
    bool APVPeakMode = true;   // to be included in header configuration file
    if (APVPeakMode == true) {
        std::ifstream inputFile(PulseFileName);

        if (!inputFile.is_open()) {
            std::cout << "Error opening the file!" << std::endl;
        }
        std::string line;
        std::string res_find = "resolution:";
        float res;
        std::string s;

        while (getline(inputFile, line)) {
            if ((!line.empty()) && (line.substr(0, 1) != "#")) {
                std::stringstream ss(line);
                if (line.find(res_find) != std::string::npos) {
                    res = stof(line.substr(line.find(res_find) + res_find.length()));   // implement check for
                                                                                        // resolution
                } else {
                    std::string value;
                    while (getline(ss, value, ' ')) {
                        PulseValues.push_back(stod(value));
                    }
                }
            }
        }
        const auto max_value = max_element(PulseValues.begin(), PulseValues.end());
        if (abs(*max_value - 1) > numeric_limits<double>::epsilon()) {
            throw invalid_argument("Maximum value of pulse shape not 1.");
        }

        unsigned int pulset0Idx = std::distance(PulseValues.begin(), max_value);
    }
    // get the vector of beginning to max of the pulse
    // time response not included!
}

TVector3 ChargeDivision::getLocal(Int_t detID, TVector3 point)
{
    TVector3 local_point; 
    // Calculate the detector id as per the geofile, where strips are disrespected
    // int strip = (detID) % 1024;                // actual strip ID
    int geofile_detID = detID; // the det id number needed to read the geometry
    int station = geofile_detID >> 17;
    int plane = (geofile_detID >> 16) % 2;
    int row = (geofile_detID >> 13) % 8;
    int column = (geofile_detID >> 11) % 4;
    int sensor = geofile_detID;
    int sensor_module = advsnd::target::columns * row + 1 + column;

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
    Double_t x = point.X();
    Double_t y = point.Y();
    Double_t z = point.Z();
    double global_pos[3] = {x, y, z};
    double local_pos[3];
    nav->MasterToLocal(global_pos, local_pos);

    local_point = local_pos; 
    return local_point; 
}



TVector3 ChargeDivision::DriftDir(TVector3 EntryPoint, TVector3 ExitPoint, float length)
{
    TVector3 DriftDirUnit = (EntryPoint - ExitPoint).Unit();
    TVector3 DriftMul = DriftDirUnit * length;
    TVector3 DriftPos = EntryPoint - (DriftDirUnit * length);
    return DriftPos;
}

EnergyFluctUnit ChargeDivision::Divide(Int_t detID, const std::vector<AdvTargetPoint*>& V)
{

    for (int i = 0; i < V.size(); i++) {

        // lorentz angle - check how it changes with radiation damaged sensors
        // need to develop logic that can calculate the lorentz angle for different temperatures, depletion voltage and
        // bias voltage !!! value below taken from "Determination of the Lorentz Angle in Microstrip Silicon Detectors
        // with Cosmic Muons" but for magnetic field of 4T!!
        //  float lorentz_angle = -4.5;
        // not required for target

        // mass, momentum, energy loss in MeV, segment length in mm
        std::vector<Double_t> fluctEnergy;
        std::vector<TVector3> driftPos;
        std::vector<TVector3> glob_driftPos;
        Int_t pdgcode = V[i]->PdgCode();


        if (TDatabasePDG::Instance()->GetParticle(pdgcode)) {
            ParticleMass = (TDatabasePDG::Instance()->GetParticle(V[i]->PdgCode())->Mass()) * 1000;
            ParticleCharge = TDatabasePDG::Instance()->GetParticle(V[i]->PdgCode())->Charge();
        } else {
            cout << "Couldn't find particle " << pdgcode << endl;
        };

        // GET THE LOCAL POSITION
        TVector3 local_entry_point = getLocal(V[i] -> GetDetectorID(), V[i]->GetEntryPoint());
        TVector3 local_exit_point = getLocal(V[i] -> GetDetectorID(), V[i]->GetExitPoint());


        int strip = (detID) % 1024;

        double len = (local_entry_point - local_exit_point).Mag();

        if (fabs(ParticleMass) < 1e-6 || ParticleCharge == 0) {
            NumberofSegments = 1;
        } else {
            NumberofSegments = 1
                               + (ChargeDivisionsperStrip
                                  * abs((local_entry_point.X() - local_exit_point.X())
                                        / StripPitch));   // check if this is the correct way
        }

        segLen = (len / NumberofSegments) * 10;   // in mm

        SiG4UniversalFluctuation sig4fluct{};

        Double_t Etotal = V[i]->GetEnergyLoss() * 1000;
        Double_t Emean = Etotal / NumberofSegments;

        // check which coordinate to consider
        Double_t momentum = sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2)) * 1000;

        if (NumberofSegments > 1) {
            for (Int_t j = 0; j < NumberofSegments; j++) {

                // sig4fluct.SampleFluctuations(ParticleMass, ParticleCharge, Emean, momentum, segLen);
                fluctEnergy.push_back(
                    sig4fluct.SampleFluctuations(ParticleMass, ParticleCharge, Emean, momentum, segLen));
                driftPos.push_back(DriftDir(local_entry_point, local_exit_point, (segLen * j) / 10));
                glob_driftPos.push_back(DriftDir(V[i]->GetEntryPoint(), V[i]->GetExitPoint(), (segLen * j) / 10));
            }
        } else {
            fluctEnergy.push_back(V[i]->GetEnergyLoss() * 1000);
            driftPos.push_back(DriftDir(local_entry_point, local_exit_point, (segLen) / 10));
            glob_driftPos.push_back(DriftDir(V[i]->GetEntryPoint(), V[i]->GetExitPoint(), (segLen) / 10));
        }

        double sume = 0;
        for (int n = 0; n < size(fluctEnergy); n++) {
            sume = sume + fluctEnergy[n];
        }
        if (sume > 0.) {
            float rescale_ratio = Etotal / sume;
            for (int m = 0; m < size(fluctEnergy); m++) {
                fluctEnergy[m] = (fluctEnergy[m] * rescale_ratio) / 1000;
            }
        }
        EnergyFluctUnit ELossVector(fluctEnergy, segLen / 10, driftPos, glob_driftPos);

        //         //to check if local position is working 
        
        
        return ELossVector;
    }
}
