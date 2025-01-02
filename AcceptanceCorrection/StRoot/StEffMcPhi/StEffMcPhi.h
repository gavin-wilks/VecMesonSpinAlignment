#ifndef StEffMcPhi_h
#define StEffMcPhi_h
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "StRoot/StEffMcPhi/StEffHistManger.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include <string>

class TFile;
class StEffHistManger;
class TChain;

class StEffMcPhi : public TObject
{
  public:
    StEffMcPhi(const char*, int Energy, long StartEvent, long StopEvent, int PID, int mode, int etamode, int inputpt, int startpt, int stoppt);
    ~StEffMcPhi();

    void SetInPutFile(const std::string inputfile);
    void SetOutPutFile(const std::string outputfile);
    void SetStartEvent(long StartEvent);
    void SetStopEvent(long StopEvent);

    void Init();
    void Make();
    void Finish();

  private:
    std::string mInPutFile;
    std::string mOutPutFile;
    long mStartEvent;
    long mStopEvent;
    TFile *mFile_InPut;
    TFile *mFile_OutPut;

    int energy;
    int pid;
    int mMode; 
    int mEtaMode;
    int mInputPt;
    int mStartPt;
    int mStopPt;

    //double mStartInputPt[3] = {1.0,2.0,3.0}; 
    //double mStopInputPt[3]  = {3.0,4.0,5.0}; 

    TChain *mInPut;
    static int mInPut_flag;
    const char* mList;
 
    StEffHistManger *mEffHistManger;

    float mCentrality;
    float mMcPt;
    float mMcP;
    float mMcEta;
    float mMcY;
    float mMcPhi;
    float mMcCos;
    float mMcCosTheta;
    float mMcKpEta;
    float mMcKmEta;
    float mMcCosRP;
    float mMcBeta;
    float mMcBetaP;

  ClassDef(StEffMcPhi,1)
};

#endif
