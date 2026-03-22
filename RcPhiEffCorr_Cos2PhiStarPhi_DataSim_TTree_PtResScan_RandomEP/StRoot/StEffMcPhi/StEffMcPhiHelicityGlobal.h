#ifndef StEffMcPhiHelicityGlobal_h
#define StEffMcPhiHelicityGlobal_h
#include "StMessMgr.h"
#include <string>
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TF2.h"

using namespace std;

class TNtuple;
class TTree;
class TFile;
class StEffHistMangerHelicityGlobal;
class StEffCut;
class TF1;
class TF2;
class TH3F;
class TH2F;
class TH1F;
class TProfile;

class StEffMcPhiHelicityGlobal
{
  public:
    StEffMcPhiHelicityGlobal(int Energy, long StartEvent, long StopEvent, int PID, int Year, int Cut, int inputpt, int startpt, int stoppt, const char* setting, int etamode, int order, int iter, int study, int method, float Rerho1n1, float Rho00, float ReTerms, float ImTerms, float Imrho1n1, float rhohelicity, float, float);
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

    void SetBinCos(int bin){ mBinCos = bin; };
    void SetBinPhi(int bin){ mBinPhi = bin; };
    void SetTofFlag(int flag){ mTofFlag = flag; };
    void SetSingleKaonFlag(int flag){ mSingleKaonFlag = flag; };

    void SetIter(int flag){ mIter = flag; };
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
    TFile *mInPutFile_Res;
  
    //double mRes[9];
    TF1 *f_res[0];
    double mChi[9];
    TF1 *f_pDel[9];

    //double mRes1[9];
    TF1 *f_res1[9];
    double mChi1[9];
    TF1 *f_pDel1[9];

    TF1* f_y;
    TF1* f_mRhoPt[vmsa::pt_rebin];
    TF1* f_mRhoPt_Helicity[vmsa::pt_rebin];
    TF1* f_mRhoCent[vmsa::pt_rebin_cent][9];
    TF1* f_mRhoY[vmsa::pt_rebin_y][vmsa::cent_rebin_total][vmsa::y_total];
    TF1* f_mV2[9]; 
    TF1* pythiaflat[19];

    TFile *ToFFile;
    TH3F* ToFHist[2];

    TFile *VzFile;
    TH1F* h_mVz; 
    TH1F* h_mVzOut; 
    TProfile *h_mCos2PhiPsi;
    TProfile *h_mCosPhiPsi;
    TProfile *h_mCos2PhiPsi1;
    TProfile *h_mCosPhiPsi1;
    TProfile *h_mV1[2];
    TProfile *h_mV2[2];
    TH2F *h_mPsiPsiRandom;

    //double mChi[2][9]; // order of EP, centrality

    int mTofFlag;

    int energy;
    int pid;
    int mInputPt;
    int mStartPt;
    int mStopPt;
    int mEtaMode;
    int mMode;
    int mOrder;
    int mMethod;

    int mIter;

    static int mInput_flag;

    int mSingleKaonFlag = 1;

    StEffCut *mEffCut;
    StEffHistMangerHelicityGlobal *mEffHistManger;

    TTree *mTTree; // event header
    Int_t mCentrality;

    Double_t mPsi;
    Double_t mPsi1;
    Double_t mPsi2;

    Double_t mPhiPx;
    Double_t mPhiPy;
    Double_t mPhiPz;
    Double_t mPhiE;

    Double_t mKpPx;
    Double_t mKpPy;
    Double_t mKpPz;
    Double_t mKpE;

    Double_t mKmPx;
    Double_t mKmPy;
    Double_t mKmPz;
    Double_t mKmE;

    float rerho1n1;
    float imrho1n1;
    float real;
    float imag;
    float mrho00helicity;
   
    float rhoiter[10][5];
    float mCenters[10];
    float mWidths[10];

    float mv2;

    float mMax = 0.0;
    float mMaxData[6];
    float mMaxHelicity[6];

    int mBinCos = 20;
    int mBinPhi = 20;
    int mStudy = 0;

    ClassDef(StEffMcPhiHelicityGlobal,1)
};

#endif
