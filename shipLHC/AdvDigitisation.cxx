#include "AdvDigitisation.h"

#include "AdvTargetPoint.h"
#include "ChargeDivision.h"
#include "EnergyFluctUnit.h"
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
    EnergyFluctUnit EnergyLossVector = chargedivision.Divide(detID, V);

    EFluct = EnergyLossVector.getEfluct();
    segLen = EnergyLossVector.getsegLen();
    DriftPos = EnergyLossVector.getDriftPos();
    glob_DriftPos = EnergyLossVector.getglobDriftPos();

    std::vector<int> temp_efluctsize;
    std::vector<std::vector<Double_t>> temp_efluct;
    std::vector<float> temp_seglen;
    std::vector<std::vector<TVector3>> temp_driftpos;
    std::vector<std::vector<TVector3>> temp_globdriftpos;

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

    ofstream myfile;
    myfile.open("digi.txt", std::ios_base::app);
    for (int i = 0; i < temp_efluctsize.size(); i++) {
        for (int j = 0; j < temp_efluct[i].size(); j++) {
            myfile << temp_efluctsize[i] << "\t" << temp_seglen[i] << "\t" << temp_efluct[i][j] << "\t" << V[i]->GetEnergyLoss() << "\t" << V[i]->PdgCode() << "\t" << sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2))<< "\t" << temp_driftpos[i][j].X() << "\t" << temp_driftpos[i][j].Y() << "\t" << temp_driftpos[i][j].Z() <<  "\t" << temp_globdriftpos[i][j].X() <<  "\t" << temp_globdriftpos[i][j].Y() <<  "\t" << temp_globdriftpos[i][j].Z() <<  endl;
        }
    }
    myfile.close();



    // if (gSystem->AccessPathName("digi.root")) {
    //     std::cout << "test does not exist" << std::endl;
    //     TFile *file = TFile::Open("digi.root", "NEW");
    //     TTree *tree = new TTree("digi", "digi");

    //     int efluctsize = 0;
    //     Double_t efluct = 0;
    //     float seglen = 0;
    //     Double_t dpX = 0;
    //     Double_t dpY = 0;
    //     Double_t dpZ = 0;
    //     Double_t eloss = 0;
    //     int num = 0;
    //     Int_t pdg = 0;
    //     Double_t mom = 0;

    //     tree->Branch("num", &num, "num/I");
    //     tree->Branch("efluctsize", &efluctsize, "efluctsize/I");
    //     tree->Branch("efluct", &efluct, "efluct/D");
    //     tree->Branch("seglen", &seglen, "seglen/F");
    //     tree->Branch("eloss", &eloss, "eloss/D");
    //     tree->Branch("pdgcode", &pdg, "pdgcode/I");
    //     tree->Branch("mom", &mom, "mom/D");
    //     tree->Branch("dpX", &dpX, "dpX/D");
    //     tree->Branch("dpY", &dpX, "dpY/D");
    //     tree->Branch("dpZ", &dpZ, "dpZ/D");

    //     for (int i = 0; i < temp_efluctsize.size(); i++) {
    //         for (int j = 0; j < temp_efluct[i].size(); j++) {
    //             num = i;
    //             efluctsize = temp_efluctsize[i];
    //             seglen = temp_seglen[i];
    //             efluct = temp_efluct[i][j];
    //             eloss = V[i]->GetEnergyLoss();
    //             pdg = V[i]->PdgCode();
    //             mom = sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2));
    //             dpX = temp_driftpos[i][j].X();
    //             dpY = temp_driftpos[i][j].Y();
    //             dpZ = temp_driftpos[i][j].Z();
    //             tree->Fill();
    //         }
    //     }
    //     // tree->Fill();
    //     tree->Write();
    //     file->Close();

    // } else {
    //     std::cout << "test exists" << std::endl;
    //     TFile fcurrent("digi.root", "UPDATE");
    //     TTree *input_tree = (TTree *)fcurrent.Get("digi");

    //     int efluctsize;
    //     Double_t efluct;
    //     float seglen;
    //     Double_t dpX, dpY, dpZ, eloss, mom;
    //     int num;
    //     Int_t pdg;

    //     input_tree->SetBranchAddress("num", &num);
    //     input_tree->SetBranchAddress("efluctsize", &efluctsize);
    //     input_tree->SetBranchAddress("efluct", &efluct);
    //     input_tree->SetBranchAddress("seglen", &seglen);
    //     input_tree->SetBranchAddress("eloss", &eloss);
    //     input_tree->SetBranchAddress("pdgcode", &pdg);
    //     input_tree->SetBranchAddress("mom", &mom);
    //     input_tree->SetBranchAddress("dpX", &dpX);
    //     input_tree->SetBranchAddress("dpY", &dpY);
    //     input_tree->SetBranchAddress("dpZ", &dpZ);

    //     for (int i = 0; i < temp_efluctsize.size(); i++) {
    //         for (int j = 0; j < temp_efluct[i].size(); j++) {
    //             num = i;
    //             efluctsize = temp_efluctsize[i];
    //             seglen = temp_seglen[i];
    //             efluct = temp_efluct[i][j];
    //             eloss = V[i]->GetEnergyLoss();
    //             pdg = V[i]->PdgCode();
    //             mom = sqrt(pow(V[i]->GetPx(), 2) + pow(V[i]->GetPy(), 2) + pow(V[i]->GetPz(), 2));
    //             dpX = temp_driftpos[i][j].X();
    //             dpY = temp_driftpos[i][j].Y();
    //             dpZ = temp_driftpos[i][j].Z();

    //             input_tree->Fill();
    //         }
    //     }
    //     input_tree->Write("", TObject::kOverwrite);
    //     fcurrent.Close();
    // }
}
