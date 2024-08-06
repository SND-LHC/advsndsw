//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4UniversalFluctuation
//
// Author:        V. Ivanchenko for Laszlo Urban
//
// Creation date: 03.01.2002
//
// Modifications:
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// #include <external/clhep/include/CLHEP/Random/RandFlat.h>
// #include <CLHEP/Random/RandGaussQ.h>
// #include <CLHEP/Random/RandPoissonQ.h>

#include "ShipUnit.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"

#include <SiG4UniversalFluctuation.h>
#include <TDatabasePDG.h>
#include <iostream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiG4UniversalFluctuation::SiG4UniversalFluctuation() {rndmarray = new Double_t[sizearray];}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void SiG4UniversalFluctuation::InitialiseMe(Int_t pdgcode)
// {
//     particleMass = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();
//     q = TDatabasePDG::Instance()->GetParticle(pdgcode)->Charge();

//     m_Inv_particleMass = 1.0 / particleMass;
//     m_massrate = ShipUnit::electron_mass_c2 * m_Inv_particleMass;
//     chargeSquare = q * q;
// }

SiG4UniversalFluctuation::~SiG4UniversalFluctuation()
{
  delete [] rndmarray;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t SiG4UniversalFluctuation::SampleFluctuations(Double_t particleMass, Double_t q, Double_t averageLoss, Double_t mom, Double_t length)
{
    // Calculate actual loss from the mean loss.
    // The model used to get the fluctuations is essentially the same
    // as in Glandz in Geant3 (Cern program library W5013, phys332).
    // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

    // shortcut for very small loss or from a step nearly equal to the range
    // (out of validity of the model)


    m_Inv_particleMass = 1.0 / particleMass;
    m_massrate = ShipUnit::electron_mass_c2 * m_Inv_particleMass;
    chargeSquare = q * q;
    

    Double_t energyLoss = 0;

    minLoss = 10. * ShipUnit::eV;
    if (averageLoss < minLoss) {
        energyLoss = averageLoss;
        return energyLoss;
    }

    meanLoss = averageLoss;

    Double_t tkin = (mom * mom) / (particleMass * particleMass) + 1.0;
    const Double_t gam = tkin * m_Inv_particleMass + 1.0;
    const Double_t gam2 = gam * gam;
    const Double_t beta = 1.0 - 1.0 / gam2;
    const Double_t beta2 = beta * beta;

    Double_t loss(0.), siga(0.);

    // Gaussian regime
    // for heavy particles only and conditions
    // for Gauusian fluct. has been changed

    tmax = 2. * ShipUnit::electron_mass_c2 * beta2 * gam2 / (1. + m_massrate * (2. * gam + m_massrate));

    if (particleMass > ShipUnit::electron_mass_c2 && meanLoss >= minNumberInteractionsBohr * tcut
        && tmax <= 2. * tcut) {

        siga =
            std::sqrt((tmax / beta2 - 0.5 * tcut) * ShipUnit::twopi_mc2_rcl2 * length * chargeSquare * electronDensity);

        const Double_t sn = meanLoss / siga;
        // cout << "siga : " << siga << " --- SN : " << sn << " --- length " << length << " --- tcut " << tcut << " --- tmax " << tmax << " --- beta2 " << beta2 << " --- m_massrate " << m_massrate << " --- gam2 " << gam2 << " --- mass " << particleMass << endl; 


        // thick target case
        if (sn >= 2.0) {

            const Double_t twomeanLoss = meanLoss + meanLoss;
            do {

                gRandom->SetSeed(0);
                TRandom* rngSaved = static_cast<TRandom*>(gRandom->Clone());
                energyLoss = gRandom->Gaus(meanLoss, siga);
                // 	// Loop checking, 03-Aug-2015, Vladimir Ivanchenko
            } while (0.0 > loss || twomeanLoss < loss);
        } 
        else {
            const Double_t neff = sn * sn;
            // check the gamma function
            // TF1* f1 = new TF1("f1", "TMath::Gamma(neff,1)");
            TF1 *f1 = new TF1("Gamma(x)","ROOT::Math::tgamma(x)",neff,1.0);
            double r = f1->GetRandom();
            gRandom->SetSeed(0);
            TRandom* rngSaved = static_cast<TRandom*>(gRandom->Clone());
            energyLoss = meanLoss * r / neff;
 
        }
        return energyLoss;
        
    }
    if (tcut <= e0) {
        return meanLoss;
    }

    const Double_t scaling = std::min(1. + 0.5 * ShipUnit::keV / tcut, 1.50);
    meanLoss /= scaling;

    w2 = (tcut > ipotFluct) ? TMath::Log(2. * ShipUnit::electron_mass_c2 * beta2 * gam2) - beta2 : 0.0;
    //return 0.33;
    return SampleGlandz() * scaling;
}

Double_t SiG4UniversalFluctuation::SampleGlandz()
{
    Double_t a1(0.0), a3(0.0);
    Double_t loss = 0.0;
    Double_t e1 = ipotFluct;
    Double_t rate = rateFluct;

    if (tcut > e1) {
        a1 = meanLoss * (1. - rate) / e1;
        if (a1 < a0) {
            const Double_t fwnow = 0.1 + (fw - 0.1) * std::sqrt(a1 / a0);
            a1 /= fwnow;
            e1 *= fwnow;
        } else {
            a1 /= fw;
            e1 *= fw;
        }
    }
    
    const Double_t w1 = tcut / e0;
    a3 = rate * meanLoss * (tcut - e0) / (e0 * tcut * TMath::Log(w1));
    if (a1 <= 0.) {
        a3 /= rate;
    }

    //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont
    Double_t emean = 0.;
    Double_t sig2e = 0.;

    // AddExcitation(a1, e1, emean, loss, sig2e);

    


    // excitation of type 1
    if (a1 > 0.0) {
        AddExcitation(a1, e1, emean, loss, sig2e);
    }

   
    if (sig2e > 0.0) {
        SampleGauss(emean, sig2e, loss);
    }
     

    // ionisation
    if (a3 > 0.) {
        emean = 0.;
        sig2e = 0.;
        Double_t p3 = a3;
        Double_t alfa = 1.;
        if (a3 > nmaxCont) {
            alfa = w1 * (nmaxCont + a3) / (w1 * nmaxCont + a3);
            const Double_t alfa1 = alfa * TMath::Log(alfa) / (alfa - 1.);
            const Double_t namean = a3 * w1 * (alfa - 1.) / ((w1 - 1.) * alfa);
            emean += namean * e0 * alfa1;
            sig2e += e0 * e0 * namean * (alfa - alfa1 * alfa1);
            p3 = a3 - namean;
        }


        const Double_t w3 = alfa * e0;
        if (tcut > w3) {
            gRandom->SetSeed(0);
            TRandom* rndm = static_cast<TRandom*>(gRandom->Clone());
            const Double_t w = (tcut - w3) / tcut;
            const Int_t nnb = (Int_t)gRandom->Poisson(p3);
            if (nnb > 0) {
                if (nnb > sizearray) {
                    sizearray = nnb;
                    delete[] rndmarray;
                    rndmarray = new Double_t[nnb];
                }
                // cout << " nnb : " << nnb << endl ;
                // // rndmEngineF->flatArray(nnb, rndmarray);
                // for (int m = 0; m < sizeof(rndmarray); ++m){cout << " entry " << m << " " << rndmarray[m] << endl; }
                // cout << " nnb : " << nnb << " --- w : " << w << " --- size : " << sizeof(rndmarray) << " --- w3 : " << w3 << endl; 

                for (Int_t k = 0; k < nnb; ++k) {
                    loss += w3 / (1. - w * rndmarray[k]);
                }
            }
        }

        if (sig2e > 0.0) {
            SampleGauss(emean, sig2e, loss);
        }
    }
    return loss;
}
