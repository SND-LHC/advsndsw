#ifndef FLUKAnuGENERATOR_H
#define FLUKAnuGENERATOR_H 1

#include "FairGenerator.h"
#include "TROOT.h"
#include "TTree.h"   // for TTree

class FairPrimaryGenerator;

class NeutronGenerator_FLUKA : public FairGenerator
{
  public:
    /** default constructor **/
    NeutronGenerator_FLUKA();

    /** destructor **/
    virtual ~NeutronGenerator_FLUKA();

    /** public method ReadEvent **/
    Bool_t ReadEvent(FairPrimaryGenerator*);
    Bool_t Init(const char*, int) { return true; };   //!
    Bool_t Init(const char*) { return true; };        //!
    void SetZ(Double_t X) { SND_Z = X; };

  private:
  protected:
    Double_t SND_Z;
    int fNevents;
    int fn;
    ClassDef(NeutronGenerator_FLUKA, 1);
};

#endif /* !FLUKAnuGENERATOR_H */
