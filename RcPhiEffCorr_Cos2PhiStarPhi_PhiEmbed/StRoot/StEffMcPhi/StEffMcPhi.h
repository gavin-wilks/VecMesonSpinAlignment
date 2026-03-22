#ifndef StEffMcPhi_h
#define StEffMcPhi_h
#include "StMessMgr.h"
#include <string>
#include <map>
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TF2.h"
#include "TLorentzVector.h"

using namespace std;

class TNtuple;
class TFile;
class StEffHistManger;
class StEffCut;
class TF1;
class TF2;
class TH1F;

typedef std::map<std::string,TF1*> TF1Map;
typedef std::map<std::string,TH1F*> TH1FMap;

class StEffMcPhi
{
  public:
    StEffMcPhi(int Energy, long StartEvent, long StopEvent, int PID, int Year, int Cut, int inputpt, int startpt, int stoppt, const char* setting, int etamode, int order);
    ~StEffMcPhi();

    void SetInPutFile(const string inputfile);
    void SetOutPutFile(const string outputfile);
    void SetStartEvent(long StartEvent);
    void SetStopEvent(long StopEvent);

    void setSigmay(std::string sigmay) { mSigmay = sigmay; }
    void setV2(float v2) { mv2 = v2; }

    void setRho00(float rho) { rho00 = rho; }
    void setReRho1n1(float rho) { rerho1n1 = rho; }
    void setImRho1n1(float rho) { imrho1n1 = rho; }
    void setReal(float rho) { real = rho; }
    void setImag(float rho) { imag = rho; }

    void setHRho00(float rho) { hrho00 = rho; }
    void setHReRho1n1(float rho) { hrerho1n1 = rho; }
    void setHImRho1n1(float rho) { himrho1n1 = rho; }
    void setHReal(float rho) { hreal = rho; }
    void setHImag(float rho) { himag = rho; }

    void findFuncPID(TLorentzVector const& lKaon, int icharge, int &EtaBin, int &PhiBin);
    bool pass_m2_PID(int icharge, int icent, TLorentzVector const& lKaon, int EtaBin, int PhiBin);
    bool pass_nsig_PID(int icharge, int icent, TLorentzVector const& lKaon, int EtaBin, int PhiBin);
    void readm2PID();
    void readnsigPID();

    void readEfficiency(int energy);
    void findHist(TLorentzVector* lPhi, int iParticleIndex, int& EtaBin, int& PhiBin);
    bool tpcReconstructed(int iParticleIndex, int cent, TLorentzVector* lPhi);
  
    void Init();
    void InitMap();
    void Make();
    void Finish();

  private:
    string mSigmay;
    string mInPutFile;
    string mOutPutFile;
    long mStartEvent;
    long mStopEvent;
    TFile *mFile_InPut;
    TFile *mFile_OutPut;
    TFile *mInPutFile_Res1;
    TFile *mInPutFile_Res2;

    TF1* f_y;
    TF1* f_rhoy;
    TF1* f_mRhoY;
    TF1* f_pDel1;
    TF1* f_pDel2;
    TF1* f_mRhoPt[vmsa::pt_rebin];
    TF1* f_mRhoPt_Helicity[vmsa::pt_rebin];
    TF2* f_mRhoPt_2D[vmsa::pt_rebin];
    TF2* f_mRhoPt_Helicity_2D[vmsa::pt_rebin];
    TF2* f_mRhoY_Helicity_2D[10];
    TF1* f_mRhoCent[vmsa::pt_rebin_cent][9];
    //TF1* f_mRhoY[vmsa::pt_rebin_y][vmsa::cent_rebin_total][vmsa::y_total];
    TF1* f_mV2[9]; 

    TF1Map f_nsig_PID;
    TF1Map f_m2_PID;

    TH1FMap h_EffPhi;
    TH1F *h_FrameEta[2];
    TH1F *h_FramePhi[2];
  
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
    StEffHistManger *mEffHistManger;

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
    float mRcTpc;

    float rho00;
    float rerho1n1;
    float imrho1n1;
    float real;
    float imag;

    float hrho00;
    float hrerho1n1;
    float himrho1n1;
    float hreal;
    float himag;
   
    float mv2;

    float mMax[6];
    float mMax2D[6];
    float mMaxHelicity[6];
    float mMaxHelicity2D[6];
    float mMaxHelicity2DY[10];

    float mMaxData[6];

    ClassDef(StEffMcPhi,1)
};

#endif
