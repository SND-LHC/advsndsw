#ifndef PNDGeGENERATOR_H
#define PNDGeGENERATOR_H 1

#include "FairGenerator.h"
#include "FairLogger.h"   // for FairLogger, MESSAGE_ORIGIN
#include "TF1.h"          // for TF1
#include "TH1.h"          // for TH1
#include "TH2.h"          // for TH2
#include "TROOT.h"
#include "TTree.h"   // for TTree
#include "TVector3.h"
#include "vector"

class FairPrimaryGenerator;

class GenieGenerator : public FairGenerator
{
  public:
    /** default constructor **/
    GenieGenerator();

    /** destructor **/
    virtual ~GenieGenerator();

    /** public method ReadEvent **/
    Bool_t OldReadEvent(FairPrimaryGenerator*);
    Bool_t ReadEvent(FairPrimaryGenerator*);
    virtual Bool_t Init(const char*, int);   //!
    virtual Bool_t Init(const char*);        //!
    Int_t GetNevents();
    void NuOnly() { fNuOnly = true; }
    void SetPositions(Double_t zTa, Double_t zS = -3400., Double_t zE = 2650.)
    {
        ztarget = zTa;
        startZ = zS;
        endZ = zE;
    }
    void AddBox(TVector3 dVec, TVector3 box);
    Double_t MeanMaterialBudget(const Double_t* start, const Double_t* end, Double_t* mparam);
    void SetDeltaE_Matching_FLUKAGenie(Double_t DeltaE) { fDeltaE_GenieFLUKA_nu = DeltaE; }
    void SetGenerationOption(Int_t GenOption) { fGenOption = GenOption; }
    void SetCrossingAngle(Double_t crossingangle) { fcrossingangle = crossingangle; }

  private:
    std::vector<double> Rotate(Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz);
    Int_t ExtractEvent_Ekin(Double_t Ekin, Double_t DeltaE);

  private:
  protected:
    Double_t FLUKA_x_cos, FLUKA_y_cos, FLUKA_x, FLUKA_y;
    Double_t fDeltaE_GenieFLUKA_nu;
    Double_t Yvessel, Xvessel, Lvessel, ztarget, startZ, endZ;
    Double_t Ev, pxv, pyv, pzv, El, pxl, pyl, pzl, vtxx, vtxy, vtxz, vtxt;
    Double_t fcrossingangle;   // crossing angle of the beam protons in LHC
    Bool_t cc, nuel;
    Int_t nf, neu;
    Int_t fGenOption;
    int fNevents;
    int fn;
    bool fFirst, fNuOnly;
    Double_t fznu0, fznu11, fXnu11, fYnu11;
    Double_t fEntrDz_inner, fEntrDz_outer, fEntrZ_inner, fEntrZ_outer, fEntrA, fEntrB, fL1z, fScintDz;

    FairLogger* fLogger;                              //!   don't make it persistent, magic ROOT command
    Double_t Ef[500], pxf[500], pyf[500], pzf[500];   //!
    Int_t pdgf[500];                                  //!
    TFile* fInputFile;                                //!
    TTree* fTree;                                     //!
    TTree* fFLUKANuTree;                              //!
    std::vector<TVector3> dVecs;                      //!
    std::vector<TVector3> boxs;                       //!
    TH1D* pxhist[3000];                               //!
    TH1D* pyslice[3000][500];                         //!
    TClonesArray* ancstr;                             //!

    ClassDef(GenieGenerator, 2);
};

#endif /* !PNDGeGENERATOR_H */
