#include "ChargeDivision.h"

#include "AdvTargetPoint.h"
#include "SiG4UniversalFluctuation.h"

#include <TDatabasePDG.h>
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
                        PulseValues.push_back(
                            stod(value));   // implement check to see if the max value of the pulse shape is 1
                    }
                }
            }
        }
    }
    // get the vector of beginning to max of the pulse
    // check if pulse correction is needed
}

std::vector<Double_t> ChargeDivision::Divide(Int_t detID, const std::vector<AdvTargetPoint*>& V)
{

    std::vector<Double_t> fluctEnergy;
    for (int i = 0; i < V.size(); i++) {

        // lorentz angle - check how it changes with radiation damaged sensors
        // need to develop logic that can calculate the lorentz angle for different temperatures, depletion voltage and
        // bias voltage !!! value below taken from "Determination of the Lorentz Angle in Microstrip Silicon Detectors
        // with Cosmic Muons" but for magnetic field of 4T!!
        //  float lorentz_angle = -4.5;
        // not required for target

        // mass, momentum, energy loss in MeV, segment length in mm

        Int_t pdgcode = V[i]->PdgCode();
        if (TDatabasePDG::Instance()->GetParticle(pdgcode)) {
            ParticleMass = (TDatabasePDG::Instance()->GetParticle(V[i]->PdgCode())->Mass()) * 1000;
            ParticleCharge = TDatabasePDG::Instance()->GetParticle(V[i]->PdgCode())->Charge();
        } else {
            cout << "Couldn't find particle " << pdgcode << endl;
        };

        double a = V[i]->GetEntryPoint().X() - V[i]->GetExitPoint().X();
        double b = V[i]->GetEntryPoint().Z() - V[i]->GetExitPoint().Z();
        double c = sqrt(pow(a, 2) + pow(b, 2));

        if (fabs(ParticleMass) < 1e-6 || ParticleCharge == 0) {
            NumberofSegments = 1;
        } else {
            NumberofSegments = ChargeDivisionsperStrip * abs(c / b);
        }

        segLen = (c / NumberofSegments) * 10;   // check why some are seglen are too small

        SiG4UniversalFluctuation sig4fluct{};

        Double_t Etotal = V[i]->GetEnergyLoss() * 1000;
        Double_t Emean = Etotal / NumberofSegments;

        // check which coordinate to consider
        Double_t momentum = V[i]->GetPx() * 1000;

        if (NumberofSegments > 1) {
            for (Int_t j = 0; j < NumberofSegments; j++) {

                sig4fluct.SampleFluctuations(ParticleMass, ParticleCharge, Emean, momentum, segLen);
                fluctEnergy.push_back(
                    sig4fluct.SampleFluctuations(ParticleMass, ParticleCharge, Emean, momentum, segLen));
            }
        } else {
            fluctEnergy.push_back(V[i]->GetEnergyLoss() * 1000);
        }

        double sume = 0;
        // Write to the file
        for (int n = 0; n < size(fluctEnergy); n++) {
            sume = sume + fluctEnergy[n];
        }
        if (sume > 0.) {
            float rescale_ratio = Etotal / sume;
            for (int m = 0; m < size(fluctEnergy); m++) {
                fluctEnergy[m] = fluctEnergy[m] * rescale_ratio;
            }
        }

        // Close the file
    }
    return fluctEnergy;
}
