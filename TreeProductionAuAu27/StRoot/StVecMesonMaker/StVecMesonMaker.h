#ifndef StVecMesonMaker_h
#define StVecMesonMaker_h

#include "StMaker.h"
#include "TString.h"
#include <iostream>
#include <fstream>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
//class StRunIdEventsDb;
class StCombPID;
class StVecMesonCut;
class StVecMesonProManger;
class StVecMesonCorrection;
class StVecMesonHistoManger;
class StVecMesonTree;
class StUtility;

class StVecMesonMaker : public StMaker {
  public:
    StVecMesonMaker(const char *name, StPicoDstMaker *picoMaker, const char *jobCounter, const Int_t Mode, const Int_t Energy, const Int_t Flag_ME, const Int_t Flag_PID);
    virtual ~StVecMesonMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent *mPicoEvent;
    static StRefMultCorr *mRefMultCorr;
    StUtility *mUtility;
    //StRunIdEventsDb *mRunIdEventsDb;
    StVecMesonCut *mVecMesonCut;
    StVecMesonProManger *mVecMesonProManger;
    StVecMesonCorrection *mVecMesonCorrection;
    StVecMesonHistoManger *mVecMesonHistoManger;
    StVecMesonTree *mVecMesonTree;
    
    Int_t mMode;
    Int_t mEnergy;
    Int_t mFlag_ME;
    Int_t mFlag_PID;


    std::ofstream mFile_OutPut;
 
    std::string mName_PID[3] = {"Phi","Rho","KStar"};

    TString mInPut_Corr_ReCenter;

    TString mOutPut_ReCenterPar;
    TString mOutPut_ShiftPar;
    TString mOutPut_Resolution;
    TString mOutPut_PID;

    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_Resolution;
    TFile *mFile_PID;

    Int_t mUsedTrackCounter;

    ClassDef(StVecMesonMaker, 1)
};

#endif
