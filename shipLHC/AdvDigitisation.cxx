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

    ChargeDivision chargedivision{};
    std::vector<EnergyFluctUnit> EnergyLossVector = chargedivision.Divide(detID, V);

    ChargeDrift chargedrift{};
    std::vector<SurfaceSignal> DiffusionSignal = chargedrift.Drift(EnergyLossVector);

    InducedCharge inducedcharge{};
    std::vector<AdvSignal> ResponseSignal = inducedcharge.IntegrateCharge(DiffusionSignal);

    // std::vector<Int_t> Strips = ResponseSignal.getStrips(); 
    // std::vector<Double_t> IntegratedSignal = ResponseSignal.getIntegratedSignal();
    // std::vector<std::vector<Double_t>> PulseResponse = ResponseSignal.getPulseResponse();

    // EFluct = EnergyLossVector.getEfluct();
    // segLen = EnergyLossVector.getsegLen();
    // DriftPos = EnergyLossVector.getDriftPos();
    // glob_DriftPos = EnergyLossVector.getglobDriftPos();

    // std::vector<int> temp_efluctsize;
    // std::vector<std::vector<Double_t>> temp_efluct;
    // std::vector<float> temp_seglen;
    // std::vector<std::vector<TVector3>> temp_driftpos;
    // std::vector<std::vector<TVector3>> temp_globdriftpos;

    // std::vector<std::vector<TVector3>> temp_diffpos;
    // std::vector<std::vector<Double_t>> temp_diffarea;
       

    // Double_t xx = DriftPos[0].X();

    // if (EFluct.empty()) {
    //     EFluctSize = 0;
    // } else {
    //     EFluctSize = size(EFluct);
    // }

    // temp_efluctsize.push_back(EFluctSize);
    // temp_efluct.push_back(EFluct);
    // temp_seglen.push_back(segLen);
    // temp_driftpos.push_back(DriftPos);
    // temp_globdriftpos.push_back(glob_DriftPos);

    // temp_diffarea.push_back(DiffusionSignal.getDiffusionArea());
    // temp_diffpos.push_back(DiffusionSignal.getSurfacePos());

    ofstream rawdatafile;
    rawdatafile.open("rawdata.txt", std::ios_base::app);

    ofstream chargedivisionfile;
    chargedivisionfile.open("chargedivision.txt", std::ios_base::app);

    ofstream chargedriftfile;
    chargedriftfile.open("chargedrift.txt", std::ios_base::app);

    ofstream inducedchargefile;
    inducedchargefile.open("inducedcharge.txt", std::ios_base::app);

    TVector3 local_entry;
    TVector3 local_exit;

    for (int i = 0; i < V.size(); i ++)
    {
        local_entry = getLocal(V[i] -> GetDetectorID(), V[i]->GetEntryPoint());
        local_exit = getLocal(V[i] -> GetDetectorID(), V[i]->GetExitPoint());
        rawdatafile << V.size() << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->PdgCode() << "\t" << sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2)) << "\t" << V[i]->GetPx() << "\t" << V[i]->GetPy() << "\t" << V[i]->GetPz() << "\t" << (V[i]->GetEntryPoint()).X() << "\t" << (V[i]->GetEntryPoint()).Y() << "\t" << (V[i]->GetEntryPoint()).Z() << "\t" << local_entry.X() << "\t" << local_entry.Y() << "\t" << local_entry.Z() << "\t" << (V[i]->GetExitPoint()).X() << "\t" << (V[i]->GetExitPoint()).Y() << "\t" << (V[i]->GetExitPoint()).Z() << "\t" << local_exit.X() << "\t" << local_exit.Y() << "\t" << local_exit.Z() << "\t" << V[i]->GetLength() << "\t" << V[i]->GetEventID() << "\t" << V[i]->GetTrackID() << "\t" << V[i]->GetTime() << "\t" << V[i]->GetDetectorID() << endl; 

        std::vector<Double_t> EFluct;
        int EFluctSize;
        float segLen;
        std::vector<TVector3> DriftPos;
        std::vector<TVector3> glob_DriftPos;

        EFluct = EnergyLossVector[i].getEfluct();
        EFluctSize = EFluct.size();
        segLen = EnergyLossVector[i].getsegLen();
        DriftPos = EnergyLossVector[i].getDriftPos();
        glob_DriftPos = EnergyLossVector[i].getglobDriftPos();
        
        std::vector<TVector3> DiffPos;
        std::vector<Double_t> DiffArea;

        DiffArea = DiffusionSignal[i].getDiffusionArea();
        DiffPos = DiffusionSignal[i].getSurfacePos();

        std::vector<Int_t> Strips; 
        std::vector<Double_t> IntegratedSignal;
        std::vector<std::vector<Double_t>> PulseResponse;

        Strips = ResponseSignal[i].getStrips(); 
        IntegratedSignal = ResponseSignal[i].getIntegratedSignal();
        PulseResponse = ResponseSignal[i].getPulseResponse();

        
        
        for (int j = 0; j < EFluct.size(); j++)
        {
            chargedivisionfile << V.size() << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->GetEventID() << "\t" << V[i]->GetTrackID() << "\t" << V[i]->GetTime() << "\t" << V[i]->GetDetectorID() << "\t" << EFluctSize << "\t" << EFluct[j] << "\t" << segLen << "\t" << DriftPos[j].X() << "\t" << DriftPos[j].Y() << "\t" << DriftPos[j].Z() << "\t" << glob_DriftPos[j].X() << "\t" << glob_DriftPos[j].Y() << "\t" << glob_DriftPos[j].Z() << endl;  

            chargedriftfile << V.size() << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->GetEventID() << "\t" << V[i]->GetTrackID() << "\t" << V[i]->GetTime() << "\t" << V[i]->GetDetectorID() << "\t" << DiffArea[j] << "\t" << DiffPos[j].X() << "\t" << DiffPos[j].Y() << "\t" << DiffPos[j].Z() << endl; 
        }

        for (int l = 0; l < Strips.size(); l++)
        {
            for (int m = 0; m < PulseResponse[l].size(); m++)
            {
            inducedchargefile << V.size() << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->GetEventID() << "\t" << V[i]->GetTrackID() << "\t" << V[i]->GetTime() << "\t" << V[i]->GetDetectorID() << "\t" << Strips[l] << "\t" << IntegratedSignal[l] << "\t" << PulseResponse[l][m] << endl; 
            }
        }
    }

    rawdatafile.close();
    chargedivisionfile.close();
    chargedriftfile.close();
    inducedchargefile.close();



    // ofstream myfile;
    // myfile.open("test.txt", std::ios_base::app);
    // for (int i = 0; i < temp_efluctsize.size(); i++) {
    //     for (int j = 0; j < temp_efluct[i].size(); j++) {
            
    //         myfile << temp_efluctsize[i] << "\t" << temp_seglen[i] << "\t" << temp_efluct[i][j] << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->PdgCode() << "\t" << sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2))<< "\t" << temp_driftpos[i][j].X() << "\t" << temp_driftpos[i][j].Y() << "\t" << temp_driftpos[i][j].Z() <<  "\t" << temp_globdriftpos[i][j].X() <<  "\t" << temp_globdriftpos[i][j].Y() <<  "\t" << temp_globdriftpos[i][j].Z() <<  "\t" << temp_diffpos[i][j].X() << "\t" << temp_diffpos[i][j].Y() << "\t" << temp_diffpos[i][j].Z() << "\t" << temp_diffarea[i][j] << endl;
    //     }
    // }
    // TVector3 local_entry = getLocal(V[0] -> GetDetectorID(), V[0]->GetEntryPoint());
    // TVector3 local_exit = getLocal(V[0] -> GetDetectorID(), V[0]->GetExitPoint());
    
    // ofstream myfile2; 
    // myfile2.open("test2.txt", std::ios_base::app);
    // for (int j = 0; j < Strips.size(); j ++)
    //     {
    //         for(int k = 0; k < PulseResponse[j].size(); k++)
    //         {
    //             myfile2 << V.size() << "\t" << temp_eloss[0] << "\t" << temp_eloss[0]/(3.61*1e-9) << "\t" << Strips[j] << "\t" << IntegratedSignal[j] << "\t" << local_entry.X() << "\t" << local_entry.Z() << "\t" << local_exit.X() << "\t" << local_exit.Z() << "\t" << PulseResponse[j][k] << endl; 
        
    //         }
    //     }
    // myfile.close();
    // myfile2.close();
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