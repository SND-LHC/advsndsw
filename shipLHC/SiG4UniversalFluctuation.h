//
// GEANT4 Class header file
//
//
// File name:   SiG4UniversalFluctuation
//
// Author: Vladimir Ivanchenko make a class for Laszlo Urban model
//
// Modified for standalone use in CMSSW. Danek K. 02/2006
//
// Class Description:
//
// Implementation of energy loss fluctuations in Silicon

// -------------------------------------------------------------------
//

#ifndef SiG4UniversalFluctuation_h
#define SiG4UniversalFluctuation_h

#include "TRandom.h"

namespace CLHEP {
class HepRandomEngine;
}

class SiG4UniversalFluctuation
{
  public:
    SiG4UniversalFluctuation();
    void InitialiseMe(Int_t pdgcode);
    Double_t SampleFluctuations(Double_t averageLoss, Double_t mom, Double_t length);
    Double_t SampleGlandz();

  protected:
    Double_t particleMass;
    Double_t q;
    Double_t m_Inv_particleMass;
    Double_t m_massrate;
    Double_t chargeSquare;

    Double_t minLoss;
    Double_t meanLoss;
    Double_t tmax;
    Double_t massrate;

    Double_t minNumberInteractionsBohr = 10.0;
    Double_t tcut = 0.120425;               // in MeV, value taken from CMS
    Double_t electronDensity = 6.797E+20;   // value taken from CMS
    Double_t e0 = 1.E-5;                    // material->GetIonisation()->GetEnergy0fluct();
    Double_t ipotFluct = 0.0001736;         // material->GetIonisation()->GetMeanExcitationEnergy();
    Double_t ipotLogFluct = -8.659;         // material->GetIonisation()->GetLogMeanExcEnergy();
    Double_t rateFluct = 0.4;               // material->GetIonisation()->GetRateionexcfluct();

    Double_t a0 = 42.0;
    Double_t fw = 4.0;
    Double_t nmaxCont = 8.0;
    Double_t w2 = 0.0;
    Int_t sizearray = 30;
    Double_t* rndmarray = nullptr;

    virtual Double_t SampleGlandz(const Double_t tcut);

    inline void AddExcitation(const Double_t ax, const Double_t ex, Double_t& eav, Double_t& eloss, Double_t& esig2);

    inline void SampleGauss(const Double_t eav, const Double_t esig2, Double_t& eloss);
};

inline void SiG4UniversalFluctuation::AddExcitation(const Double_t ax,
                                                    const Double_t ex,
                                                    Double_t& eav,
                                                    Double_t& eloss,
                                                    Double_t& esig2)
{
    if (ax > nmaxCont) {
        eav += ax * ex;
        esig2 += ax * ex * ex;
    } else {
        gRandom->SetSeed(0);
        TRandom* rngSaved = static_cast<TRandom*>(gRandom->Clone());
        const Int_t p = (Int_t)gRandom->Poisson(ax);
        if (p > 0) {
            eloss += ((p + 1) - 2. * gRandom->Uniform()) * ex;
        }
    }
}

inline void SiG4UniversalFluctuation::SampleGauss(const Double_t eav, const Double_t esig2, Double_t& eloss)
{
    Double_t x = eav;
    const Double_t sig = std::sqrt(esig2);
    gRandom->SetSeed(0);
    TRandom* rndm = static_cast<TRandom*>(gRandom->Clone());
    if (eav < 0.25 * sig) {
        x += (2. * rndm->Uniform() - 1.) * eav;
    } else {
        do {
            x = rndm->Gaus(eav, sig);
        } while (x < 0.0 || x > 2 * eav);
        // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
    }
    eloss += x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
