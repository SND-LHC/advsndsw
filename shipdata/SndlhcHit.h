#ifndef SNDLHCHIT_H
#define SNDLHCHIT_H 1

#include "Rtypes.h"   // for Double_t, Int_t, Double32_t, etc
#include "TArrayF.h"
#include "TObject.h"   //
#include "TVector3.h"

#include <unordered_map>

#ifndef __CINT__
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#endif   //__CINT__

/**
 * copied from FairRoot FairHit and simplified
 */
class SndlhcHit : public TObject
{

  public:
    /** Default constructor **/
    SndlhcHit();

    /** Constructor with detector id, number of SiPMs per side, number of sides **/
    SndlhcHit(Int_t detID, Int_t nSiPMs = 1, Int_t nSides = 0);

    /** Destructor **/
    virtual ~SndlhcHit();

    /** Accessors **/
    Int_t GetDetectorID() const { return fDetectorID; };
    Float_t GetSignal(Int_t nChannel = 0);
    Float_t GetTime(Int_t nChannel = 0);
    Int_t GetnSiPMs() const { return nSiPMs; };
    Int_t GetnSides() const { return nSides; };
    /** Modifiers **/
    void SetDigi(Float_t s, Float_t t, Int_t i = 0)
    {
        signals[i] = trunc(100 * s) / 100;
        times[i] = trunc(1000 * t) / 1000.;
    }
    void SetDetectorID(Int_t detID) { fDetectorID = detID; }
    void SetDaqID(Int_t i, Int_t k, Int_t board_id, Int_t tofpet_id, Int_t tofpet_channel)
    {
        fDaqID[i] = k * 100000 + board_id * 1000 + tofpet_id * 100 + tofpet_channel;
    }
    Int_t GetBoardID(Int_t i) { return int((fDaqID[i] % 100000) / 1000); }
    Int_t GetTofpetID(Int_t i) { return int((fDaqID[i] % 1000) / 100); }
    Int_t Getchannel(Int_t i) { return fDaqID[i] % 100; }
    Int_t GetRawHitIndex(Int_t i = 0) { return int(fDaqID[i] / 100000); }

    // to be implemented by the subdetector

    /*** Output to screen */
    virtual void Print(const Option_t* opt = "") const { ; }
    /*** Get position */
    virtual void GetPosition(TVector3 L, TVector3 R) const { ; }
    /*** Get energy */
    // virtual Float_t GetEnergy();  // causes problems, don't know why: cling::DynamicLibraryManager::loadLibrary():
    // lib/libShipData.so.0.0.0: undefined symbol: _ZN9SndlhcHit9GetEnergyEv

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& boost::serialization::base_object<TObject>(*this);
        ar & fDetectorID;
        ar & nSiPMs;
        ar & nSides;
    }

  protected:
#ifndef __CINT__   // for BOOST serialization
    friend class boost::serialization::access;
#endif                     // for BOOST serialization
    Int_t fDetectorID;     ///< Detector unique identifier
    Int_t nSiPMs;          /// number of SiPMs per side
    Int_t nSides;          /// number of sides
    Float_t signals[16];   /// SiPM signal
    Float_t times[16];     /// SiPM time
    Int_t fDaqID[16];      /// encodes rawhitindex*100000+(board_id * 1000) + (tofpet_id * 100) + tofpet_channel

    ClassDef(SndlhcHit, 2);
};

#endif
