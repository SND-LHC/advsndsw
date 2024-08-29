#include "AdvDigitisation.h"

#include "AdvTargetPoint.h"
#include "ChargeDivision.h"
#include "ChargeDrift.h"
#include "InducedCharge.h"
#include "EnergyFluctUnit.h"
#include "SurfaceSignal.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

AdvDigitisation::AdvDigitisation() {}

void AdvDigitisation::digirun(Int_t detID, const std::vector<AdvTargetPoint *> &V)
{

    std::vector<Int_t> temp_detID ; //-->need to change this
    for(int i =0; i < V.size(); i++)
    {
        temp_detID.push_back(V[i]->GetDetectorID());
    }
    ChargeDivision chargedivision{};
    EnergyFluctUnit EnergyLossVector = chargedivision.Divide(detID, V);

    ChargeDrift chargedrift{};
    //chargedrift.Drift(EnergyLossVector);
    SurfaceSignal DiffusionSignal = chargedrift.Drift(EnergyLossVector);

    InducedCharge inducedcharge{};
    inducedcharge.IntegrateCharge(temp_detID, DiffusionSignal);


    EFluct = EnergyLossVector.getEfluct();
    segLen = EnergyLossVector.getsegLen();
    DriftPos = EnergyLossVector.getDriftPos();
    glob_DriftPos = EnergyLossVector.getglobDriftPos();

    std::vector<int> temp_efluctsize;
    std::vector<std::vector<Double_t>> temp_efluct;
    std::vector<float> temp_seglen;
    std::vector<std::vector<TVector3>> temp_driftpos;
    std::vector<std::vector<TVector3>> temp_globdriftpos;

    std::vector<std::vector<TVector3>> temp_diffpos;
    std::vector<std::vector<Double_t>> temp_diffarea;
       

    Double_t xx = DriftPos[0].X();

    if (EFluct.empty()) {
        EFluctSize = 0;
    } else {
        EFluctSize = size(EFluct);
    }

    temp_efluctsize.push_back(EFluctSize);
    temp_efluct.push_back(EFluct);
    temp_seglen.push_back(segLen);
    temp_driftpos.push_back(DriftPos);
    temp_globdriftpos.push_back(glob_DriftPos);

    temp_diffarea.push_back(DiffusionSignal.getDiffusionArea());
    temp_diffpos.push_back(DiffusionSignal.getSurfacePos());

    ofstream myfile;
    myfile.open("digi.txt", std::ios_base::app);
    for (int i = 0; i < temp_efluctsize.size(); i++) {
        for (int j = 0; j < temp_efluct[i].size(); j++) {
            
            myfile << temp_efluctsize[i] << "\t" << temp_seglen[i] << "\t" << temp_efluct[i][j] << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->PdgCode() << "\t" << sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2))<< "\t" << temp_driftpos[i][j].X() << "\t" << temp_driftpos[i][j].Y() << "\t" << temp_driftpos[i][j].Z() <<  "\t" << temp_globdriftpos[i][j].X() <<  "\t" << temp_globdriftpos[i][j].Y() <<  "\t" << temp_globdriftpos[i][j].Z() <<  "\t" << temp_diffpos[i][j].X() << "\t" << temp_diffpos[i][j].Y() << "\t" << temp_diffpos[i][j].Z() << "\t" << temp_diffarea[i][j] << endl;
        }
    }
    myfile.close();
}
