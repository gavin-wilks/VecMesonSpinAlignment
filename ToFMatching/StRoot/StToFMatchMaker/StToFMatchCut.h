#ifndef StToFMatchCut_h
#define StToFMatchCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class StRefMultCorr;
class StPicoEvent;

class StToFMatchCut : public TObject
{
  public:
    StToFMatchCut(Int_t energy);
    virtual ~StToFMatchCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*, StPicoEvent*);
    bool passSigPionCut(StPicoTrack*, StPicoEvent*, Float_t);
    bool passSigProntonCut(StPicoTrack*, StPicoEvent*, Float_t);
    bool passSigKaonCut(StPicoTrack*, StPicoEvent*, Float_t);
    bool passToFMatchCut(StPicoTrack*, StPicoEvent*, StPicoDst*);
    bool isMinBias(StPicoEvent*);
    bool isPileUpEvent(int, int, double);
    double getRefMultReweight(double, int);
    int getCentrality(double);
    double getEventWeight(int, double);
    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();
    Float_t getMass2(StPicoTrack*, StPicoDst*);
    Float_t getV0Mass2(StPicoTrack*, StPicoDst*);

  private:
    static StRefMultCorr *mRefMultCorr;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;
    Int_t mEnergy;

    ClassDef(StToFMatchCut,1)
};
#endif
