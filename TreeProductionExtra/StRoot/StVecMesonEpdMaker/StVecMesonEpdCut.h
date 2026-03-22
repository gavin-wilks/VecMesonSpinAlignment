#ifndef StVecMesonEpdCut_h
#define StVecMesonEpdCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

class StVecMesonEpdCut : public TObject
{
  public:
    StVecMesonEpdCut(const int energy);
    virtual ~StVecMesonEpdCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*, StPicoEvent*);
    bool passSigPionCut(StPicoTrack*, Int_t);
    bool passSigKaonCut(StPicoTrack*, Int_t);
    bool passTrackMeson(StPicoTrack*, StPicoEvent*, Int_t);
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
   
    ClassDef(StVecMesonEpdCut,1)
};
#endif
