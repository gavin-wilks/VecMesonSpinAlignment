#ifndef StVecMesonAna_h
#define StVecMesonAna_h

#include "TObject.h"
#include "TString.h"
#include "THn.h"
//#include "../Utility/phi_data_constants_19GeV.h"

class StRefMultCorr;
class TFile;
class TChain;
class StMesonEvent;
class StMesonTrack;
class StVecMesonCorr;
class StVecMesonCut;
class StVecMesonHistoManger;
//class StRunIdEventsDb;
class StUtility;
class TTree;

class StVecMesonAna : public TObject
{
  public:
    StVecMesonAna(const Char_t *list, const Char_t *jobId, Int_t energy, Int_t X_flag, Int_t mode, Int_t etamode); // X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for rho, 2 for Kstar
    ~StVecMesonAna();

    void setInputDir(const TString inputdir);
    void setOutputfile(const TString outputfile);
    void setInPutList(const TString iInPutList);

    void Init();
    void InitSE();
    void InitME();
    void Make();
    void MakePhi();
    void Finish();

  private:
    TString mInputdir;
    TString mOutputfile;
    TString mInPutList;

    TTree *mKaonTreeSE;
    int   mCentSE; 
    float mWeightSE; 
    int   mChargeSE;
    float mPtSE;
    float mRapiditySE;
    float mEtaSE;
    float mPhiSE;
    Bool_t mHasTofInfoSE;

    TTree *mKaonTreeME;
    int   mCentME; 
    float mWeightME; 
    int   mChargeME;
    float mPtME;
    float mRapidityME;
    float mEtaME;
    float mPhiME;
    Bool_t mHasTofInfoME;

    //int   mNHitsFit;
    //int   mNHitsMax;
    //float mDEdx;
    //float mDca;

    TFile *mFile_OutPut;
    TChain *mInPutSE;
    TChain *mInPutME;
    Int_t mEnergy;
    Int_t mEtaMode;
    float mEtaCut;
    Int_t mX_flag; // 0 for Same Event, 1 for Mixed Event
    const Char_t *mList;
    Int_t mMode; // 0 for phi, 1 for rho, 2 for Kstar
    const Char_t *mJobId;
    StMesonEvent *mMeson_eventSE;
    StMesonEvent *mMeson_eventME;
    StMesonTrack *mMeson_track;
    StVecMesonCorr *mVecMesonCorr;
    StVecMesonCut *mVecMesonCut;
    StVecMesonHistoManger *mVecMesonHistoManger;
    //StRunIdEventsDb *mRunIdEventsDb;
    StUtility *mUtility;  

    THnF *ptyetaphi[2]; 

    static StRefMultCorr *mRefMultCorr;
    static Int_t mInPut_flag;
    static char* VM_EVENT_TREE_SE;
    static char* VM_EVENT_TREE_ME;
    static char* VM_EVENT_BRANCH;

  ClassDef(StVecMesonAna,1)
};
#endif
