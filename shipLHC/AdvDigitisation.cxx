#include "AdvDigitisation.h"
#include "AdvTargetPoint.h"
#include "ChargeDivision.h"
#include "ChargeDrift.h"
#include "InducedCharge.h"
#include "EnergyFluctUnit.h"
#include "SurfaceSignal.h"
#include "AdvSignal.h"
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

    std::vector<Double_t> temp_eloss ; //-->need to change this
    for(int i =0; i < V.size(); i++)
    {
        temp_eloss.push_back(V[i]->GetEnergyLoss());
    }
    ChargeDivision chargedivision{};
    EnergyFluctUnit EnergyLossVector = chargedivision.Divide(detID, V);

    ChargeDrift chargedrift{};
    //chargedrift.Drift(EnergyLossVector);
    SurfaceSignal DiffusionSignal = chargedrift.Drift(EnergyLossVector);

    InducedCharge inducedcharge{};
    AdvSignal ResponseSignal = inducedcharge.IntegrateCharge(DiffusionSignal);

    std::vector<Int_t> Strips = ResponseSignal.getStrips(); 
    std::vector<Double_t> IntegratedSignal = ResponseSignal.getIntegratedSignal();

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
    myfile.open("test.txt", std::ios_base::app);
    for (int i = 0; i < temp_efluctsize.size(); i++) {
        for (int j = 0; j < temp_efluct[i].size(); j++) {
            
            myfile << temp_efluctsize[i] << "\t" << temp_seglen[i] << "\t" << temp_efluct[i][j] << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->PdgCode() << "\t" << sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2))<< "\t" << temp_driftpos[i][j].X() << "\t" << temp_driftpos[i][j].Y() << "\t" << temp_driftpos[i][j].Z() <<  "\t" << temp_globdriftpos[i][j].X() <<  "\t" << temp_globdriftpos[i][j].Y() <<  "\t" << temp_globdriftpos[i][j].Z() <<  "\t" << temp_diffpos[i][j].X() << "\t" << temp_diffpos[i][j].Y() << "\t" << temp_diffpos[i][j].Z() << "\t" << temp_diffarea[i][j] << endl;
        }
    }
    TVector3 local_entry = getLocal(V[0] -> GetDetectorID(), V[0]->GetEntryPoint());
    TVector3 local_exit = getLocal(V[0] -> GetDetectorID(), V[0]->GetExitPoint());
    
    ofstream myfile2; 
    myfile2.open("test2.txt", std::ios_base::app);
    for (int j = 0; j < Strips.size(); j ++)
        {
            myfile2 << V.size() << "\t" << temp_eloss[0] << "\t" << temp_eloss[0]/(3.61*1e-9) << "\t" << Strips[j] << "\t" << IntegratedSignal[j] << "\t" << local_entry.X() << "\t" << local_entry.Z() << "\t" << local_exit.X() << "\t" << local_exit.Z() << endl; 
        }
    myfile.close();
    myfile2.close();
}

TVector3 AdvDigitisation::getLocal(Int_t detID, TVector3 point)
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