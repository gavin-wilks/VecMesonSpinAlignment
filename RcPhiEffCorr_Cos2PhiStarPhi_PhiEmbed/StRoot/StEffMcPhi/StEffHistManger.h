#ifndef StEffHistManger_h
#define StEffHistManger_h
#include "TObject.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/phi_data_constants_19GeV.h"
#include "StRoot/Utility/type.h"
#include <map>
#include <string>

class TH1D;
class TH1F;
class TH2D;
class TH2F;
class TH3D;

typedef std::map<TString,TH2F*> TH2FMap;
typedef std::map<TString,TH1F*> TH1FMap1D;

class StEffHistManger : public TObject
{
  public:
    StEffHistManger(int energy, int pid, int mode, int startpt, int stoppt);
    virtual ~StEffHistManger();
    void InitHist();
    void InitPhiHist();
    void InitKaonHist();
    void FillPhiHistMc(int cent, float pt, float phi, float y, float phistar, float cos, float cosH, float kpt, float ky);
    void FillPhiHistRc(int cent, float pt, float phi, float y, float phistar, float cos, float cosH, float kpt, float ky);
    void FillKaonHistMc(int cent, float phi, float pt, float y, float kpt, float ky, float keta, float kphi, float phistar, float cos, float cosH);
    void FillKaonHistRc(int cent, float phi, float pt, float y, float kpt, float ky, float keta, float kphi, float phistar, float cos, float cosH);
    void FillKaonDeltaHistMc(int cent, float phi, float eta, float pt, float y, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar);
    void FillKaonDeltaHistRc(int cent, float phi, float eta, float pt, float y, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar);
    void FillHistMc(int,float,float,float,float,float,float,float,float,float,float,float,float,float);
    void FillHistRc(int,float,float,float,float,float,float,float,float,float,float,float,float,float,int);
    float AngleShift(float);
    //void CalEfficiency();
    //void CalEffPtEtaPhi();
    void CalEffCosThetaStar();
    TH1D* CalEffError(TH1D*,TH1D*,std::string);
    TH2D* CalEffError(TH2D*,TH2D*,std::string);
    void WriteHist();
    void WritePhiHist();
    void WriteKaonHist();

  private:
    TH2FMap h_mPhi_MC;
    TH2FMap h_mPhi_RC;
    TH1FMap1D h_mPhi1D_MC;
    TH1FMap1D h_mPhi1D_RC;
    TH2FMap h_mKaon_MC;
    TH2FMap h_mKaon_RC;
    TH2F *h_mKaon_MC_Array[2][10];
    TH2F *h_mKaon_RC_Array[10][2][10];

    TH3D *h_mMcTracks[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3D *h_mRcTracks[10][10];

    TH1D *h_mMcEffPt[10]; // pt distritbution as a function of centrality
    TH1D *h_mRcEffPt[10];
    TH1D *h_mEffPt[10];

    TH1D *h_mMcEffEta[10]; // pt distritbution as a function of centrality and eta
    TH1D *h_mRcEffEta[10];
    TH1D *h_mEffEta[10];

    TH1D *h_mMcEffPhi[10]; // pt distritbution as a function of centrality and phi
    TH1D *h_mRcEffPhi[10];
    TH1D *h_mEffPhi[10];

    TH1DMap h_mMcEffPEP;
    TH1DMap h_mRcEffPEP;
    TH1DMap h_mEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi

    TH1D *h_mMcEffCos[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCos[10][10][10];
    TH1D *h_mEffCos[10][10][10];

    TH2D *h_mMcEffCosCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH2D *h_mRcEffCosCosH[10][10][10];
    TH2D *h_mEffCosCosH[10][10][10];

    TH1D *h_mMcEffCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosH[10][10][10];
    TH1D *h_mEffCosH[10][10][10];

    TH1D *h_mMcEffPhiS[10][10]; // efficiency vs phistar-phi as a function of centrality and pt
    TH1D *h_mRcEffPhiS[10][10][10];
    TH1D *h_mEffPhiS[10][10][10];

    TH1D *h_mMcEffCosY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosY[10][10][vmsa::y_total];
    TH1D *h_mEffCosY[10][10][vmsa::y_total];

    TH1D *h_mMcEffPhiSY[10][10][vmsa::y_total]; // efficiency vs phistar-phi as a function of centrality and pt
    TH1D *h_mRcEffPhiSY[10][10][vmsa::y_total];
    TH1D *h_mEffPhiSY[10][10][vmsa::y_total];

    TH2D *h_mMcEffCosEP[10][10]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mRcEffCosEP[10][10][10];
    TH2D *h_mEffCosEP[10][10][10];

    TH2D *h_mMcEffCosPhiPrime[10][10]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mRcEffCosPhiPrime[10][10][10];
    TH2D   *h_mEffCosPhiPrime[10][10][10];

    TH2D *h_mMcEffCosHPhiPrime[10][10]; // helicity efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mRcEffCosHPhiPrime[10][10][10];
    TH2D   *h_mEffCosHPhiPrime[10][10][10];

    TH2D *h_mMcEffCosEPY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mRcEffCosEPY[10][10][vmsa::y_total];
    TH2D *h_mEffCosEPY[10][10][vmsa::y_total];

    TH2D *h_mMcEffPtY[10]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mRcEffPtY[10][10];
    TH2D *h_mEffPtY[10][10];

    TH2D *h_mSpectraRatio;
    TH2D *h_mKaonSpectraRatio[4];
    TH2D *h_mCosCosHRatio[4];

    //MODE 4 Rapidity dependence
    TH1D *h_mYMcEffCos[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mYRcEffCos[10][10][10];
    TH1D *h_mYEffCos[10][10][10];

    TH2D *h_mYMcEffCosCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH2D *h_mYRcEffCosCosH[10][10][10];
    TH2D *h_mYEffCosCosH[10][10][10];

    TH1D *h_mYMcEffCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mYRcEffCosH[10][10][10];
    TH1D *h_mYEffCosH[10][10][10];

    TH1D *h_mYMcEffPhiS[10][10]; // efficiency vs phistar-phi as a function of centrality and pt
    TH1D *h_mYRcEffPhiS[10][10][10];
    TH1D *h_mYEffPhiS[10][10][10];

    TH2D *h_mYMcEffCosEP[10][10]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mYRcEffCosEP[10][10][10];
    TH2D *h_mYEffCosEP[10][10][10];

    TH2D *h_mYMcEffCosPhiPrime[10][10]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mYRcEffCosPhiPrime[10][10][10];
    TH2D *h_mYEffCosPhiPrime[10][10][10];

    TH2D *h_mYMcEffCosHPhiPrime[10][10]; // helicity efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mYRcEffCosHPhiPrime[10][10][10];
    TH2D *h_mYEffCosHPhiPrime[10][10][10];


 
    int mEnergy;
    int mMode;
    int mStartPt;
    int mStopPt;

    int flag_eff;
    int flag_eff_PtEtaPhi;
    int flag_eff_Cos;
    
    int mpt_first;
    int mpt_last;

    float mpt_low[10];
    float mpt_up[10];
 
  ClassDef(StEffHistManger,1)
};

#endif
