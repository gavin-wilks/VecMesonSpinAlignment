#ifndef StToFMatchCut_h
#define StToFMatchCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class StRefMultCorr;

class StToFMatchCut : public TObject
{
  public:
    StToFMatchCut(Int_t energy);
    virtual ~StToFMatchCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passSigPionCut(StPicoTrack*, Float_t);
    bool passSigProntonCut(StPicoTrack*, Float_t);
    bool passSigKaonCut(StPicoTrack*, Float_t);
    bool passToFMatchCut(StPicoTrack*);
    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();
    Float_t getMass2(StPicoTrack*);
    Float_t getV0Mass2(StPicoTrack*);

  private:
    static StRefMultCorr *mRefMultCorr;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;
    Int_t mEnergy;

    ClassDef(StToFMatchCut,1)
};
#endif
