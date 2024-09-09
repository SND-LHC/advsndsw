#pragma once

#include "TObject.h"
#include "TVector3.h"
#include "Track.h"

class sndRecoTrack : public TObject
{
  public:
    sndRecoTrack() { ; }
    /* Constructor from genfit::Track object */
    sndRecoTrack(genfit::Track*);
    ~sndRecoTrack() { ; }

    TVector3 getStart() { return start; }
    TVector3 getStop() { return stop; }
    TVector3 getTrackMom() { return fTrackMom; }
    std::vector<TVector3> getTrackPoints() { return fTrackPoints; }
    std::vector<int> getRawMeasDetIDs() { return fRawMeasDetID; }
    std::vector<std::vector<float>> getRawMeasTimes() { return fRawMeasTimes; }
    int getTrackType() { return fTrackType; }
    bool getTrackFlag() { return fFlag; }
    float getChi2() { return chi2; }
    int getNdf() { return Ndf; }
    float getChi2Ndf() { return chi2 / (Ndf + 1E-10); }

    void setRawMeasTimes(std::vector<std::vector<float>> l) { fRawMeasTimes = l; }
    void setTrackType(int t) { fTrackType = t; }

    /* extrapolateToPlaneAtZ assumes a plane in the physics coordinate
       system is perpendicular to z-axis, which is not true for
       TI18 geometry since detector planes are tilted! */
    TVector3 extrapolateToPlaneAtZ(float z);
    std::pair<int, float> TrackDirection();
    std::pair<float, float> Velocity();
    std::tuple<float, float, float> trackDir();
    std::vector<float> getCorrTimes();
    float getSlopeXZ() { return fTrackMom.X() / fTrackMom.Z(); }
    float getSlopeYZ() { return fTrackMom.Y() / fTrackMom.Z(); }
    float getAngleXZ() { return atan(fTrackMom.X() / fTrackMom.Z()); }
    float getAngleYZ() { return atan(fTrackMom.Y() / fTrackMom.Z()); }

  private:
    std::vector<TVector3> fTrackPoints;
    std::vector<int> fRawMeasDetID;
    std::vector<std::vector<float>> fRawMeasTimes;
    TVector3 fTrackMom;
    int fTrackType;   // to be done as enum
    TVector3 start, stop;
    // selected fit status items
    bool fFlag;   // True if track fit converged
    float chi2;
    int Ndf;

    ClassDef(sndRecoTrack, 4);
};
