#ifndef StVecMesonAnaEpd_h
#define StVecMesonAnaEpd_h

#include "TObject.h"
#include "TString.h"

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

class StVecMesonAnaEpd : public TObject
{
  public:
    StVecMesonAnaEpd(const Char_t *list, const Char_t *jobId, Int_t energy, Int_t X_flag, Int_t mode, Int_t etamode); // X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for rho, 2 for Kstar
    ~StVecMesonAnaEpd();

    void setInputDir(const TString inputdir);
    void setOutputfile(const TString outputfile);
    void setInPutList(const TString iInPutList);
    void setStopEvent(const Long64_t StopEvent);
    void setStartEvent(const Long64_t StartEvent);

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
    Long64_t mStopEvent;
    Long64_t mStartEvent;

    Long64_t mStart_Event;
    Long64_t mStop_Event;

    TFile *mFile_OutPut;
    TChain *mInPut;
    Int_t mEnergy;
    Int_t mX_flag; // 0 for Same Event, 1 for Mixed Event
    const Char_t *mList;
    Int_t mMode; // 0 for phi, 1 for rho, 2 for Kstar
    Int_t mEtaMode;
    float mEtaCut;
    const Char_t *mJobId;
    StMesonEvent *mMeson_event;
    StMesonTrack *mMeson_track;
    StVecMesonCorr *mVecMesonCorr;
    StVecMesonCut *mVecMesonCut;
    StVecMesonHistoManger *mVecMesonHistoManger;
    //StRunIdEventsDb *mRunIdEventsDb;
    StUtility *mUtility;  
 
    static StRefMultCorr *mRefMultCorr;
    static Int_t mInPut_flag;
    static char* VM_EVENT_TREE;
    static char* VM_EVENT_BRANCH;

  ClassDef(StVecMesonAnaEpd,1)
};
#endif
