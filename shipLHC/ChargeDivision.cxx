#include "ChargeDivision.h"

#include "AdvTargetPoint.h"
#include "SiG4UniversalFluctuation.h"

#include <TDatabasePDG.h>
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
        // std::vector<double> pulsevalues;
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

std::vector<float> ChargeDivision::Divide(Int_t detID, const std::vector<AdvTargetPoint*>& V)
{
    for (int i = 0; i < V.size(); i++) {
        Int_t pdgcode = V[i]->PdgCode();
        double ParticleMass = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();
        double ParticleCharge = TDatabasePDG::Instance()->GetParticle(pdgcode)->Charge();

        // lorentz angle - check how it changes with radiation damaged sensors
        // need to develop logic that can calculate the lorentz angle for different temperatures, depletion voltage and
        // bias voltage !!! value below taken from "Determination of the Lorentz Angle in Microstrip Silicon Detectors
        // with Cosmic Muons" but for magnetic field of 4T!!
        //  float lorentz_angle = -4.5;
        // not required for target

        // check if this is correct
        // to be put in config files
        // calculating number of strips WITHOUT lorentz angle
        float NumberofStrips;
        // find correct pitch value
        float StripPitch = 80e-6;
        float StripWidth = 0.25 * StripPitch;
        Int_t ChargeDivisionsperStrip = 10;
        Int_t NumberofSegments = 0;
        if (ParticleCharge == 0) {
            NumberofSegments = 1;
        } else {
            NumberofStrips = ceil(V[i]->GetEntryPoint().X() - V[i]->GetExitPoint().X()) / (StripWidth + StripPitch);
            NumberofSegments = ChargeDivisionsperStrip * NumberofStrips;
        }
        // cout << NumberofSegments << endl;

        Double_t seglen = V[i]->GetLength() / NumberofSegments;

        SiG4UniversalFluctuation sig4fluct{};
        sig4fluct.InitialiseMe(pdgcode);

        Double_t El = V[i]->GetEnergyLoss() / NumberofSegments;

        // check which coordinate to consider
        Double_t momentum = V[i]->GetPx();
        // Double_t fluctEnergy[NumberofSegments];
        std::vector<float> fluctEnergy;

        for (Int_t i = 0; i < NumberofSegments; i++) {

            fluctEnergy.push_back(sig4fluct.SampleFluctuations(El, momentum, seglen));
        }
        return fluctEnergy;
    }
}
