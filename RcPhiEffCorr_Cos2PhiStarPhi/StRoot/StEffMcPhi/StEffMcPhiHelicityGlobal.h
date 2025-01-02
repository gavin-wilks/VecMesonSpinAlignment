#ifndef StEffMcPhiHelicityGlobal_h
#define StEffMcPhiHelicityGlobal_h
#include "StMessMgr.h"
#include <string>
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TF2.h"

using namespace std;

class TNtuple;
class TFile;
class StEffHistMangerHelicityGlobal;
class StEffCut;
class TF1;
class TF2;
class TH3F;

class StEffMcPhiHelicityGlobal
{
  public:
    StEffMcPhiHelicityGlobal(int Energy, long StartEvent, long StopEvent, int PID, int Year, int Cut, int inputpt, int startpt, int stoppt, const char* setting, int etamode, int order, float Rerho1n1, float Rho00, float ReTerms, float ImTerms, float Imrho1n1, float rhohelicity, float, float);
    ~StEffMcPhiHelicityGlobal();

    void SetInPutFile(const string inputfile);
    void SetOutPutFile(const string outputfile);
    void SetStartEvent(long StartEvent);
    void SetStopEvent(long StopEvent);

    void setSigmay(std::string sigmay) { mSigmay = sigmay; }
    //void setRho00(float rho) { mrho00 = rho; }
    void setV2(float v2) { mv2 = v2; }
  
    void Init();
    void InitMap();
    void Make();
    void Finish();

  private:
    string mSigmay;
    float mrho00;
    string mInPutFile;
    string mOutPutFile;
    long mStartEvent;
    long mStopEvent;
    TFile *mFile_InPut;
    TFile *mFile_OutPut;
    TFile *mInPutFile_Res1;
    TFile *mInPutFile_Res2;

    TF1* f_y;
    TF1* f_pDel1;
    TF1* f_pDel2;
    TF1* f_mRhoPt[vmsa::pt_rebin];
    TF1* f_mRhoPt_Helicity[vmsa::pt_rebin];
    TF2* f_mRhoPt_2D[5][5];
    TF1* f_mRhoCent[vmsa::pt_rebin_cent][9];
    TF1* f_mRhoY[vmsa::pt_rebin_y][vmsa::cent_rebin_total][vmsa::y_total];
    TF1* f_mV2[9]; 

    TFile *ToFFile;
    TH3F* ToFHist[2];

    double mChi[2][9]; // order of EP, centrality

    int energy;
    int pid;
    int mInputPt;
    int mStartPt;
    int mStopPt;
    int mEtaMode;
    int mMode;
    int mOrder;

    static int mInput_flag;

    StEffCut *mEffCut;
    StEffHistMangerHelicityGlobal *mEffHistManger;

    TNtuple *mNtuple; // event header
    float mCentrality;
    float mPsi;
    float mPsi1;
    float mPsi2;
    float mMcPt;
    float mMcP;
    float mMcEta;
    float mMcY;
    float mMcPhi;
    float mMcInvMass;
    float mMcPid;

    float mKpMcPt;
    float mKpMcEta;
    float mKpMcY;
    float mKpMcPhi;
    float mKpMcM;
    float mKpMcPid;

    float mKpRcPt;
    float mKpRcEta;
    float mKpRcY;
    float mKpRcPhi;
    float mKpRcM;
    float mKpRcTpc;
    float mKpRcTof;

    float mKmMcPt;
    float mKmMcEta;
    float mKmMcY;
    float mKmMcPhi;
    float mKmMcM;
    float mKmMcPid;

    float mKmRcPt;
    float mKmRcEta;
    float mKmRcY;
    float mKmRcPhi;
    float mKmRcM;
    float mKmRcTpc;
    float mKmRcTof;

    float mRcPt;
    float mRcP;
    float mRcEta;
    float mRcY;
    float mRcPhi;
    float mRcInvMass;

    float rerho1n1;
    float imrho1n1;
    float real;
    float imag;
    float mrho00helicity;
   
    float mv2;

    float mMax = 0.0;
    float mMaxData[6];
    float mMaxHelicity[6];

    ClassDef(StEffMcPhiHelicityGlobal,1)
};

#endif
