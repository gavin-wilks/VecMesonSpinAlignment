#ifndef StRunIdMatching_h
#define StRunIdMatching_h

#include "TObject.h"
#include "TString.h"
#include <map>

class StRefMultCorr;
class TFile;
class TChain;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StVecMesonCorr;
class StRunIdEventsDb;

class StRunIdMatching : public TObject
{
  public:
    StRunIdMatching(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode); // X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for Kstar, 2 for K0S
    ~StRunIdMatching();

    void setInputDir(const TString inputdir);
    void setOutputfile(const TString outputfile);
    void setInPutList(const TString iInPutList);
    void setStopEvent(const Long64_t StopEvent);
    void setStartEvent(const Long64_t StartEvent);

    void Init();
    void InitMap();
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

    std::map<std::string,Int_t> mEmbdRunId;
    std::map<std::string,Int_t> mEmbdEvtId;

    TChain *mInPut;
    Int_t mEnergy;
    Int_t mX_flag; // 0 for Same Event, 1 for Mixed Event
    Int_t mList;
    Int_t mMode; // 0 for phi, 1 for Kstar
    StAlexPhiMesonEvent *mPhiMeson_event;
    StAlexPhiMesonTrack *mPhiMeson_track;
    StVecMesonCorr *mVecMesonCorr;
    StRunIdEventsDb *mRunIdEventsDb;

    static StRefMultCorr *mRefMultCorr;
    static Int_t mInPut_flag;
    static char* VM_EVENT_TREE;
    static char* VM_EVENT_BRANCH;

  ClassDef(StRunIdMatching,1)
};
#endif
