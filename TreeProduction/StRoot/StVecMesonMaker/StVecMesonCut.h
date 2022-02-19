#ifndef StVecMesonCut_h
#define StVecMesonCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

class StVecMesonCut : public TObject
{
  public:
    StVecMesonCut(const int energy);
    virtual ~StVecMesonCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*, StPicoEvent*);
    bool passSigPionCut(StPicoTrack*, Float_t);
    bool passSigKaonCut(StPicoTrack*, Float_t);
    bool passSigProntonCut(StPicoTrack*, Float_t);
    bool passTrackMeson(StPicoTrack*, StPicoEvent*);
    bool isMinBias(StPicoEvent*);
    bool isPileUpEvent(int, int, double);
    double getRefMultReweight(double, int);
    int getCentrality(double);
    double getEventWeight(int, double);

    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();
    Float_t getBeta(StPicoTrack*, StPicoDst*);
    Float_t getPrimaryMass2(StPicoTrack*, StPicoDst*);
    Float_t getV0Mass2(StPicoTrack*, StPicoDst*);

  private:
    static StRefMultCorr *mRefMultCorr;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;
    Int_t mEnergy;
   
    ClassDef(StVecMesonCut,1)
};
#endif
