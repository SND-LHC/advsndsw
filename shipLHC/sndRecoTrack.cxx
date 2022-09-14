#include "sndRecoTrack.h"
#include "Scifi.h"
#include "MuFilter.h"
#include "ShipUnit.h"
#include "TROOT.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "StateOnPlane.h"  /* genfit */
#include "FitStatus.h"
#include "RKTrackRep.h"
#include "ConstField.h"
#include "FieldManager.h"
#include "TGeoMaterialInterface.h"
#include "MaterialEffects.h"

using namespace genfit;
using namespace std;

//tests
#include <iostream>

/* Constructor from genfit::Track object */
sndRecoTrack::sndRecoTrack(Track* track)
{
    for ( auto i = 0; i < track->getNumPoints(); i++ )
    {
         auto state = track->getFittedState(i);
         fTrackPoints.push_back(state.getPos());
         fRawMeasDetID.push_back(track->getPointWithMeasurement(i)->getRawMeasurement()->getDetId());
         if (i ==0)
         {
             fTrackMom = state.getMom();
             start = state.getPos();
         }
         if (i == track->getNumPoints()-1)
             stop = state.getPos();
   }
   FitStatus* fitStatus = track->getFitStatus();
   chi2  = fitStatus->getChi2();
   Ndf   = fitStatus->getNdf();
   fFlag = fitStatus->isFitConverged();
   // defaults
   fTrackType = 0;
   fRawMeasTimes = {};
}

pair<int, float> sndRecoTrack::TrackDirection()
{
   /* Don't use this function for now. It works ok with MC,
      not the case for data.. Work in progress. 
      Provide track direction with a probability level
      based on timing measurements of hits/clusters.
      What follows is based on a note by C.Vilela.
      Reference T0 is 1st meas in z
      1st item in the returned pair is: 
      1 = along Z axis; -1 = reverse Z direction
   */

   // Account for signal propagation in detectors
   vector<float> corr_times = getCorrTimes();

   MuFilter *MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter") );
   Scifi *ScifiDet = dynamic_cast<Scifi*> (gROOT->GetListOfGlobals()->FindObject("Scifi") );
   double resol{}, resol_0{};
   int size = fTrackPoints.size();
   TMatrixD Tres1(size-1,1), Tres2(size-1,1);
   TMatrixD Sigma(size-1, size-1);
   // Extract time meas. resolution of the reference T0 measurement
   if (fRawMeasDetID[0] >= 100000) {
         resol_0 = ScifiDet->GetConfParF("Scifi/timeResol");
   }
   else resol_0 = MuFilterDet->GetConfParF("MuFilter/timeResol");
   // Set time residual vectors and covariance matrix
   for ( int i = 0; i < size-1; i++){
       Tres1[i][0] = corr_times[i+1] - corr_times[0]
                     -(fTrackPoints[i+1]-fTrackPoints[0]).Mag()/ShipUnit::c_light;
       Tres2[i][0] = corr_times[0] - corr_times[i+1]
                     -(fTrackPoints[i+1]-fTrackPoints[0]).Mag()/ShipUnit::c_light;
       //tests
       /*cout<<i<<" deltaT "<<corr_times[i+1] - corr_times[0]<<" dL "
             << (fTrackPoints[i+1]-fTrackPoints[0]).Mag()
       <<" dL/c "<<  (fTrackPoints[i+1]-fTrackPoints[0]).Mag()/ShipUnit::c_light<<endl;*/
       for ( int j = 0; j < size-1; j++){
           if (i == j) {
              if (fRawMeasDetID[i+1] >= 100000) {
                 resol = ScifiDet->GetConfParF("Scifi/timeResol");
              }
              else resol = MuFilterDet->GetConfParF("MuFilter/timeResol");
              Sigma[i][j] = pow(resol_0,2)+pow(resol,2);
           }
           else Sigma[i][j] = pow(resol_0,2);
       }
   }
   TMatrixD SigmaInv(TMatrixD::kInverted,Sigma);
   TMatrixD Tres1_T(TMatrixD::kTransposed,Tres1);
   //tests
   //Tres1_T.Print();
   TMatrixD Tres2_T(TMatrixD::kTransposed,Tres2);
   TMatrixT<double> LambdaB1 = Tres1_T*SigmaInv*Tres1;
   TMatrixT<double> LambdaB2 = Tres2_T*SigmaInv*Tres2;
   // tests
   /*cout<<LambdaB1[0][0]<<" "<<LambdaB2[0][0]
       <<" "<<LambdaB2[0][0]-LambdaB1[0][0]
       <<" "<<TMath::Prob(LambdaB1[0][0], size-1)
       <<" "<<TMath::Prob(LambdaB2[0][0], size-1)<<endl;*/
   if (LambdaB2[0][0]-LambdaB1[0][0] >= 0)
      return make_pair(1, TMath::Prob(LambdaB1[0][0], size-1));
   else return make_pair(-1, TMath::Prob(LambdaB2[0][0], size-1));
}

pair<float, float> sndRecoTrack::Velocity()
{
   /* Extract velocity along z based on timing
      measurements as slope of linear fit of (dT,dZ)
      Reference T0 is 1st measurement in z */

   // Account for signal propagation in detectors
   vector<float> corr_times = getCorrTimes();

   MuFilter *MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter") );
   Scifi *ScifiDet = dynamic_cast<Scifi*> (gROOT->GetListOfGlobals()->FindObject("Scifi") );
   TGraph *gr = new TGraph();
   double resol{}, delta_T{};
   for (int i = 1;  i < fTrackPoints.size(); i++){
     if (fRawMeasDetID[i] >= 100000) {
         resol = ScifiDet->GetConfParF("Scifi/timeResol");
     }
     else resol = MuFilterDet->GetConfParF("MuFilter/timeResol");
     delta_T = corr_times[i]-corr_times[0];
     // Use current track point measurement only if time difference is above resolution
     if (fabs(delta_T) < resol*sqrt(2)) continue;
     gr->AddPoint(delta_T, fTrackPoints[i].Z()-fTrackPoints[0].Z());
   }
   TF1 *line = new TF1("line", "pol1");
   // A single entry in the graph - no fit
   if (gr->GetN() < 2)
     return make_pair(9999.,999.);
   else
     gr->Fit("line", "SQ");
     return make_pair(line->GetParameter(1), line->GetParError(1));
}

pair<float, float> sndRecoTrack::trackDir()
{
   /* Based on the same function in SndlhcTracking.py!
      Extract direction along z based on timing
      measurements as slope of linear fit of (dT,dZ)
      Featuring time-of-flight correction btw measurements.
      Reference T0 is 1st measurement in z */

   // Account for signal propagation in detectors
   vector<float> corr_times = getCorrTimes();

   MuFilter *MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter") );
   Scifi *ScifiDet = dynamic_cast<Scifi*> (gROOT->GetListOfGlobals()->FindObject("Scifi") );
   TGraph *gr = new TGraph();
   double resol{}, dist_Z{}, delta_T{};
   for (int i = 1;  i < fTrackPoints.size(); i++){
     if (fRawMeasDetID[i] >= 100000) {
         resol = ScifiDet->GetConfParF("Scifi/timeResol");
     }
     else resol = MuFilterDet->GetConfParF("MuFilter/timeResol");
     delta_T = corr_times[i]-corr_times[0];
     // Use current track point measurement only if time difference is above resolution
     if (fabs(delta_T) < resol*sqrt(2)) continue;
     dist_Z = fTrackPoints[i].Z()-fTrackPoints[0].Z();
     gr->AddPoint(dist_Z, delta_T - dist_Z/ShipUnit::c_light);
   }
   TF1 *line = new TF1("line", "pol1");
   // A single entry in the graph - no fit
   if (gr->GetN() < 2)
     return make_pair(9999.,999.);
   else
     gr->Fit("line", "SQ");
     return make_pair(line->GetParameter(1),
                      line->GetParameter(1)/(line->GetParError(1)+1E-13));
}

TVector3 sndRecoTrack::extrapolateToPlaneAtZ(float z)
{
   /* Set mandatory items for genfit::Extrapolate* methods
      No magnetic field and assuming no (negligible) multiple scattering */
   ConstField *bfield = new ConstField(0, 0, 0);
   FieldManager *fM = (FieldManager*)fM->getInstance();
   fM->init(bfield);
   TGeoMaterialInterface *geoMat = new TGeoMaterialInterface();
   MaterialEffects *matEff = (MaterialEffects*)matEff->getInstance();
   matEff->init(geoMat);
   matEff->setNoEffects();
   
   TVector3 NewPosition = TVector3(0., 0., z);
   /* line below assumes that plane in global coordinates
      is perpendicular to z-axis, which is not true for TI18 geometry. */
   TVector3 parallelToZ = TVector3(0., 0., 1.);
   RKTrackRep *rep = new RKTrackRep(13);
   StateOnPlane state = StateOnPlane(rep);
   // find index(!) of track point closest in z
   float zMin = 9999;
   int index;
   for ( int i =0; i < fTrackPoints.size(); i++ ){
      if ( fabs(zMin - fTrackPoints[i].Z()) < zMin ) index = i;
   }
   TVector3 Closest_pos = fTrackPoints[index];
   rep->setPosMom(state, Closest_pos, fTrackMom);
   rep->extrapolateToPlane(state, NewPosition, parallelToZ);

   return state.getPos();
}

vector<float> sndRecoTrack::getCorrTimes()
{
   /* Account for signal propagation in detectors */
   vector<float> corr_times{};
   MuFilter *MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter") );
   Scifi *ScifiDet = dynamic_cast<Scifi*> (gROOT->GetListOfGlobals()->FindObject("Scifi") );
   TVector3 A, B, X{};
   float scintVel;
   float calibTime{};
   for ( int i = 0;  i < fTrackPoints.size(); i++ ){
      if (fRawMeasDetID[i] >= 100000) {
         ScifiDet->GetSiPMPosition(fRawMeasDetID[i],A,B);
         scintVel = ScifiDet->GetConfParF("Scifi/signalSpeed");
         calibTime = ScifiDet->GetCorrectedTime(fRawMeasDetID[i], fRawMeasTimes[i], 0);
      }
      else { 
         MuFilterDet->GetPosition(fRawMeasDetID[i],A,B);
         if (floor(fRawMeasDetID[i]/10000) == 3)
            scintVel = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
         else scintVel = MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");
         // For the moment only SciFi is time calibrated
         calibTime = fRawMeasTimes[i];
      }
      // vertical detector elements
      if ( (fRawMeasDetID[i] >= 100000 && int(fRawMeasDetID[i]/100000)%10 == 1) or
           (fRawMeasDetID[i] < 100000  && floor(fRawMeasDetID[i]/10000) == 3 
                                       && fRawMeasDetID[i]%1000 > 59) ) X = B-fTrackPoints[i];
      else X = A - fTrackPoints[i];
      corr_times.push_back(fRawMeasTimes[i] - X.Mag()/scintVel);
      //cout<<"i "<<i<<" raw tdc "<<fRawMeasTimes[i]<<" corr t "<<corr_times.back()<<" "<< X.Mag()<<endl;
   }

   return corr_times;
}

