#pragma once

#include "TObject.h"
#include "TVector3.h"

class sndRecoTrack : public TObject {
 public :
  sndRecoTrack(){;}
  ~sndRecoTrack(){;}

  TVector3 getStart() {return start;}
  TVector3 getStop() {return stop;}
  std::vector<TVector3 > getTrackPoints() {return fTrackPoints;}
  std::vector<TVector3 > getTrackPointsMom() {return fTrackPointsMom;}
  std::vector<int> getRawMeasDetIDs() {return fRawMeasDetID;}
  std::vector<float> getRawMeasTimes() {return fRawMeasTimes;}
  int   getTrackType() { return fTrackType; }
  bool  getTrackFlag() { return fFlag; }
  float getChi2()      { return chi2; }
  float getChi2Ndf()   { return chi2Ndf; }

  void setStart(TVector3 s){ start = s; }
  void setStop(TVector3 s){ stop = s; }
  void setRawMeasDetIDs(std::vector<int> l) { fRawMeasDetID = l; }
  void setRawMeasTimes(std::vector<float> l) { fRawMeasTimes = l; }
  void setTrackPoints(std::vector<TVector3 > l) { fTrackPoints = l; }
  void setTrackPointsMom(std::vector<TVector3 > l) { fTrackPointsMom = l; }
  void setTrackType(int t) { fTrackType = t; }
  void setTrackFlag(bool f) { fFlag = f; }
  void setChi2(float f) { chi2 = f; }
  void setChi2Ndf(float f) { chi2Ndf = f; }

  TVector3 extrapolateToPlaneAtZ(float z);
  float TrackDirection();
  bool IsAlongZ();
  std::vector<float> getCorrTimes();

 private :
  std::vector<TVector3 > fTrackPoints;
  std::vector<TVector3 > fTrackPointsMom;
  std::vector<int> fRawMeasDetID;
  std::vector<float> fRawMeasTimes;
  int fTrackType; // to be done as enum
  TVector3 start, stop;
  // selected fit status items
  bool  fFlag; // True if track fit converged
  float chi2;
  float chi2Ndf;

  ClassDef(sndRecoTrack, 1);
};
