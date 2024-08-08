#include "ChargeDivision.h"

#include "AdvTargetPoint.h"
#include "EnergyFluctUnit.h"
#include "SiG4UniversalFluctuation.h"
#include "TVector3.h"

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
        Int_t pdgcode = V[i]->PdgCode();
        std::vector<Double_t> dx;
        std::vector<Double_t> dy;
        std::vector<Double_t> dz;

        if (TDatabasePDG::Instance()->GetParticle(pdgcode)) {
            ParticleMass = (TDatabasePDG::Instance()->GetParticle(V[i]->PdgCode())->Mass()) * 1000;
            ParticleCharge = TDatabasePDG::Instance()->GetParticle(V[i]->PdgCode())->Charge();
        } else {
            cout << "Couldn't find particle " << pdgcode << endl;
        };

        // GET THE LOCAL POSITION

        int strip = (detID) % 1024;

        double len = (V[i]->GetEntryPoint() - V[i]->GetExitPoint()).Mag();

        if (fabs(ParticleMass) < 1e-6 || ParticleCharge == 0) {
            NumberofSegments = 1;
        } else {
            NumberofSegments = 1
                               + (ChargeDivisionsperStrip
                                  * abs((V[i]->GetEntryPoint().X() - V[i]->GetExitPoint().X())
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
                driftPos.push_back(DriftDir(V[i]->GetEntryPoint(), V[i]->GetExitPoint(), (segLen * j) / 20));
                dx.push_back(driftPos[j].X());
                dy.push_back(driftPos[j].Y());
                dz.push_back(driftPos[j].Z());
            }
        } else {
            fluctEnergy.push_back(V[i]->GetEnergyLoss() * 1000);
            driftPos.push_back(DriftDir(V[i]->GetEntryPoint(), V[i]->GetExitPoint(), (segLen) / 20));
            dx.push_back(V[i]->GetExitPoint().X());
            dy.push_back(V[i]->GetExitPoint().Y());
            dz.push_back(V[i]->GetExitPoint().Z());
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
        EnergyFluctUnit ELossVector(fluctEnergy, segLen / 10, driftPos);
        return ELossVector;
    }
}
