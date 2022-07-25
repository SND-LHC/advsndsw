#include "sndRecoTrack.h"
#include "Scifi.h"
#include "MuFilter.h"
#include "ShipUnit.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "StateOnPlane.h"
#include "RKTrackRep.h"

using namespace genfit;

//tests
#include <iostream>
using namespace std;

float sndRecoTrack::TrackDirection()
{
   /* Extract direction based on timing measurements - following note by Cris
      Reference T0 is 1st meas in z */

   // Account for signal propagation in detectors
   vector<float> corr_times = getCorrTimes();
   
   MuFilter *MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter") );
   Scifi *ScifiDet = dynamic_cast<Scifi*> (gROOT->GetListOfGlobals()->FindObject("Scifi") );
   double Lambda{}, resol{};
   for ( int i = 1;  i < fTrackPoints.size(); i++ ){
      if (fRawMeasDetID[i] >= 100000) {
         resol = ScifiDet->GetConfParF("Scifi/timeResol");
      }
      else resol = MuFilterDet->GetConfParF("MuFilter/timeResol");
      Lambda+= pow(corr_times[i] - corr_times[0]- (fTrackPoints[i]-fTrackPoints[0]).Mag()/ShipUnit::c_light, 2)/pow(resol,2);
      // tests
      /*cout<<"deltaT "<<corr_times[i] - corr_times[0]<<" dL "<< (fTrackPoints[i]-fTrackPoints[0]).Mag()
      <<" dL/c "<<  (fTrackPoints[i]-fTrackPoints[0]).Mag()/ShipUnit::c_light
      << "  L[i] "<<pow(corr_times[i] - corr_times[0]- (fTrackPoints[i]-fTrackPoints[0]).Mag()/ShipUnit::c_light, 2)/pow(resol,2)<<endl;
      cout<<"accum L "<<Lambda<<endl;*/
   }
   
   return TMath::Prob(Lambda, fTrackPoints.size()-1);
}

bool sndRecoTrack::IsAlongZ()
{
   /* Extract direction based on timing measurements as
      slope of linear fit of (dT,dL)
      Reference T0 is 1st measurement in z */

   // Account for signal propagation in detectors
   vector<float> corr_times = getCorrTimes();

   TGraph *gr = new TGraph();  
   for (int i = 1;  i < fTrackPoints.size(); i++){
     gr->SetPoint(i,corr_times[i]-corr_times[0],(fTrackPoints[i]-fTrackPoints[0]).Mag());
   }
   TF1 *line = new TF1("line", "pol1");
   gr->Fit("line", "SQ");  
   if (line->GetParameter(1) >= 0) return true;
   else return false;
}

TVector3 sndRecoTrack::extrapolateToPlaneAtZ(float z)
{
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
   TVector3 Closest_mom = fTrackPointsMom[index];
   rep->setPosMom(state, Closest_pos, Closest_mom);
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
   for ( int i = 0;  i < fTrackPoints.size(); i++ ){
      if (fRawMeasDetID[i] >= 100000) {
         ScifiDet->GetSiPMPosition(fRawMeasDetID[i],A,B);
         scintVel = ScifiDet->GetConfParF("Scifi/signalSpeed");
      }
      else { 
         MuFilterDet->GetPosition(fRawMeasDetID[i],A,B);
         if (floor(fRawMeasDetID[i]/10000) == 3)
            scintVel = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
         else scintVel = MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");
      }
      // vertical detector elements
      if ( (fRawMeasDetID[i] >= 100000 && int(fRawMeasDetID[i]/100000)%10 == 1) or
           (fRawMeasDetID[i] < 100000  && floor(fRawMeasDetID[i]/10000) == 3 
                                     && fRawMeasDetID[i]%1000 > 59) ) X = B-fTrackPoints[i];
      else X = A - fTrackPoints[i];
      corr_times.push_back(fRawMeasTimes[i] - X.Mag()/scintVel);
   }

   return corr_times;
}

