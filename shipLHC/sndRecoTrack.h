#pragma once

#include "TObject.h"
#include "TVector3.h"
#include "Track.h"

class sndRecoTrack : public TObject {
 public :
  sndRecoTrack(){;}
  /* Constructor from genfit::Track object */
  sndRecoTrack(genfit::Track*);
  ~sndRecoTrack(){;}

  TVector3 getStart() {return start;}
  TVector3 getStop() {return stop;}
  TVector3 getTrackMom() {return fTrackMom;}
  std::vector<TVector3 > getTrackPoints() {return fTrackPoints;}
  std::vector<int> getRawMeasDetIDs() {return fRawMeasDetID;}
  std::vector<float> getRawMeasTimes() {return fRawMeasTimes;}
  int   getTrackType() { return fTrackType; }
  bool  getTrackFlag() { return fFlag; }
  float getChi2()      { return chi2; }
  int   getNdf()       { return Ndf; }
  float getChi2Ndf()   { return chi2/(Ndf+1E-10); }

  void setRawMeasTimes(std::vector<float> l) { fRawMeasTimes = l; }
  void setTrackType(int t) { fTrackType = t; }

  TVector3 extrapolateToPlaneAtZ(float z);
  std::pair<int, float> TrackDirection();
  std::pair<float, float> Velocity();
  std::pair<float, float> trackDir();
  std::vector<float> getCorrTimes();

 private :
  std::vector<TVector3 > fTrackPoints;
  std::vector<int> fRawMeasDetID;
  std::vector<float> fRawMeasTimes;
  TVector3 fTrackMom;
  int fTrackType; // to be done as enum
  TVector3 start, stop;
  // selected fit status items
  bool  fFlag; // True if track fit converged
  float chi2;
  int Ndf;

  ClassDef(sndRecoTrack, 3);
};
