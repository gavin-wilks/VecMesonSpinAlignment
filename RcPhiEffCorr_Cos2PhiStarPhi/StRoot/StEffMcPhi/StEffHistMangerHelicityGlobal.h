#ifndef StEffHistMangerHelicityGlobal_h
#define StEffHistMangerHelicityGlobal_h
#include "TObject.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
//#include "StRoot/Utility/phi_data_constants_19GeV.h"
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

class StEffHistMangerHelicityGlobal : public TObject
{
  public:
    StEffHistMangerHelicityGlobal(int energy, int pid, int mode, int startpt, int stoppt, int, int);
    virtual ~StEffHistMangerHelicityGlobal();
    void InitHist();
    void InitPhiHist();
    void InitKaonHist();
    void FillPhiHistMc(int cent, float pt, float phi, float y, float phistar, float cos, float cosH, float kpt, float ky);
    void FillPhiHistRc(int cent, float pt, float phi, float y, float phistar, float cos, float cosH, float kpt, float ky);
    void FillKaonHistMc(int cent, float phi, float pt, float y, float kpt, float ky, float keta, float kphi, float phistar, float cos, float cosH);
    void FillKaonHistRc(int cent, float phi, float pt, float y, float kpt, float ky, float keta, float kphi, float phistar, float cos, float cosH);
    void FillKaonDeltaHistMc(int cent, float phi, float eta, float pt, float y, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar);
    void FillKaonDeltaHistRc(int cent, float phi, float eta, float pt, float y, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar);
    void FillHistMc(int,float,float,float,float,float,float,float,float,float,float,float,float,int,int,int,int,int);
    void FillHistRc(int,float,float,float,float,float,float,float,float,float,float);
    float AngleShift(float);
    void CalEfficiency();
    void CalEffPtEtaPhi();
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

    TH3D *h_mMcTracks[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3D *h_mRcTracks[10];

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
    TH1D *h_mRcEffCos[10][10];
    TH1D *h_mEffCos[10][10];

    TH2D *h_mMcEffCosCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH2D *h_mRcEffCosCosH[10][10];
    TH2D *h_mEffCosCosH[10][10];

    TH1D *h_mMcEffCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosH[10][10];
    TH1D *h_mEffCosH[10][10];

    TH1D *h_mMcEffPhiS[10][10]; // efficiency vs phistar-phi as a function of centrality and pt
    TH1D *h_mRcEffPhiS[10][10];
    TH1D *h_mEffPhiS[10][10];

    TH1D *h_mMcEffCosY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosY[10][10][vmsa::y_total];
    TH1D *h_mEffCosY[10][10][vmsa::y_total];

    TH1D *h_mMcEffPhiSY[10][10][vmsa::y_total]; // efficiency vs phistar-phi as a function of centrality and pt
    TH1D *h_mRcEffPhiSY[10][10][vmsa::y_total];
    TH1D *h_mEffPhiSY[10][10][vmsa::y_total];

    TH2D *h_mMcEffCosEP[10][10]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mRcEffCosEP[10][10];
    TH2D *h_mEffCosEP[10][10];

    TH2D *h_mMcEffCosPhiPrime[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mMcEffCosPhiPrimePsi[11][5][5][11][20]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mRcEffCosPhiPrime[11][5];
    TH2D   *h_mEffCosPhiPrime[11][5];

    TH2D *h_mMcEffCosPhiPrimeH[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mMcEffCosPhiPrimeHPsi[11][5][5][11][20]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2D *h_mRcEffCosPhiPrimeH[11][5];
    TH2D   *h_mEffCosPhiPrimeH[11][5];

    TH3D *h3_mMcEffCosPhiPrime[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeH[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt

    TH2D *h_mMcEffCosEPY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mRcEffCosEPY[10][10][vmsa::y_total];
    TH2D *h_mEffCosEPY[10][10][vmsa::y_total];

    TH2D *h_mMcEffPtY[10]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2D *h_mRcEffPtY[10];
    TH2D *h_mEffPtY[10];

    TH2D *h_mSpectraRatio;
    TH2D *h_mKaonSpectraRatio[4];
    TH2D *h_mCosCosHRatio[4];
    TH2D *h_mV2[11];
 
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
 
  ClassDef(StEffHistMangerHelicityGlobal,1)
};

#endif
