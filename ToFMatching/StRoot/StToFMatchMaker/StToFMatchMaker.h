#ifndef StToFMatchMaker_h
#define StToFMatchMaker_h

#include "StMaker.h"
#include "TString.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StCombPID;
class StToFMatchCut;
class StToFMatchHistoManger;
class StUtility;

class StToFMatchMaker : public StMaker {
  public:
    StToFMatchMaker(const char *name, StPicoDstMaker *picoMaker, const char *jobCounter, const int Energy, const int pid);
    virtual ~StToFMatchMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent *mPicoEvent;
    static StRefMultCorr *mRefMultCorr;
    StToFMatchCut *mToFMatchCut;
    StToFMatchHistoManger *mToFMatchHistoManger;
    float AngleShift(float phi);
    
    StUtility *mUtility;

    int mEnergy;
    int mPID;

    TString mOutPut_ToFMatch;

    TFile *mFile_ToFMatch;

    Int_t mUsedTrackCounter;
    float mScaleFactor_nSigma;

    ClassDef(StToFMatchMaker, 1)
};

#endif
