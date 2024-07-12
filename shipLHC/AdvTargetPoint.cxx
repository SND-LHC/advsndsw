#include "AdvTargetPoint.h"

#include <iostream>
using std::cout;
using std::endl;

// -----   Default constructor   -------------------------------------------
AdvTargetPoint::AdvTargetPoint()
    : FairMCPoint()
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------

AdvTargetPoint::AdvTargetPoint(Int_t trackID,
                               Int_t detID,
                               TVector3 pos,
                               TVector3 mom,
                               Double_t tof,
                               Double_t length,
                               Double_t eLoss,
                               Int_t pdg_code,
                               TVector3 exit_point)
    : FairMCPoint(trackID, detID, pos, mom, tof, length, eLoss)
    , fPdgCode(pdg_code)
    , fExitX(exit_point.X())
    , fExitY(exit_point.Y())
    , fExitZ(exit_point.Z())
{
}

// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
AdvTargetPoint::~AdvTargetPoint() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void AdvTargetPoint::Print(const Option_t* opt) const
{
    cout << "-I- AdvTargetPoint: box point for track " << fTrackID << " in detector " << fDetectorID << endl;
    cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
    cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
    cout << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss " << fELoss * 1.0e06 << " keV"
         << endl;
}
// -------------------------------------------------------------------------
