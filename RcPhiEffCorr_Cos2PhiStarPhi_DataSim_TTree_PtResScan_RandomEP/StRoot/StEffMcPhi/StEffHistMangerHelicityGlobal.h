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
class TProfile;

typedef std::map<TString,TH2F*> TH2FMap;
typedef std::map<TString,TH1F*> TH1FMap1D;
typedef std::map<TString,TH2D*> TH2DMap;
typedef std::map<TString,TH1D*> TH1DMap1D;

class StEffHistMangerHelicityGlobal : public TObject
{
  public:
    StEffHistMangerHelicityGlobal(int energy, int pid, int mode, int startpt, int stoppt, int, int);
    virtual ~StEffHistMangerHelicityGlobal();
    void InitHist(int,int,int);
    void InitPhiHist();
    void InitKaonHist();
    void FillPhiHistMc(int cent, double pt, double phi, double y, double phistar, double cos, double cosH, double kpt, double ky);
    void FillPhiHistRc(int cent, double pt, double phi, double y, double phistar, double cos, double cosH, double kpt, double ky);
    void FillKaonHistMc(int cent, double phi, double pt, double y, double kpt, double ky, double keta, double kphi, double phistar, double cos, double cosH);
    void FillKaonHistRc(int cent, double phi, double pt, double y, double kpt, double ky, double keta, double kphi, double phistar, double cos, double cosH);
    void FillKaonDeltaHistMc(int cent, double phi, double eta, double pt, double y, double kppt, double kpeta, double kpphi, double kmpt, double kmeta, double kmphi, double phistar);
    void FillKaonDeltaHistRc(int cent, double phi, double eta, double pt, double y, double kppt, double kpeta, double kpphi, double kmpt, double kmeta, double kmphi, double phistar);
    void FillHistMc(int,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int);
    void FillHistCutSmearMC(int cent, double pt, double y, double cos, double phiprime, double weight, int cut, int iep, int iptres);
    void FillHistCutSmear(int cent, double pt, double y, double cos, double phiprime, double weight, int cut, int iep, int iptres);
    void FillAngleSmear(int cent, double y, double cosrc, double cosmc, double phiprimerc, double phiprimemc, double weight, int irhoinput, int irho, int cut, int iep, int iptres);
    void FillHistRandom(int cut, int cent, double pt, double y, double cosep, double cosrand, double weight);
    void FillHistCut(int,double,double,double,double,double,double,double,double,double,double,double,double,int,int,int,int,int,int);
    void FillPhiMassMC(int,double,double);
    void FillPhiMass(int,double,int,double);
    void FillPtRes(int,int,double,double,double,int,double);
    void FillHistRc(int,double,double,double,double,double,double,double,double,double,double);
    double AngleShift(double);
    void CalEfficiency();
    void CalEffPtEtaPhi();
    void CalEffCosThetaStar();
    TH1D* CalEffError(TH1D*,TH1D*,std::string);
    TH2D* CalEffError(TH2D*,TH2D*,std::string);
    void WriteHist();
    void WritePhiHist();
    void WriteKaonHist();

  private:
    TH2DMap h_mPhi_MC;
    TH2DMap h_mPhi_RC;
    TH1DMap1D h_mPhi1D_MC;
    TH1DMap1D h_mPhi1D_RC;
    TH2DMap h_mKaon_MC;
    TH2DMap h_mKaon_RC;

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

    //TH3D *h3_mMcEffCosPhiPrime[11][5][5][11][9][25]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrime[5][5][11][9][25]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeH[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeY[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeHY[11][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt

    TH3D *h3_mMcEffCosPhiPrimeYCuts[10][5][5][11][9]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeYCuts_SmearMC[2]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeYCuts_Smear[42]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mMcEffCosPhiPrimeYCuts_Smear_Cut[5][42]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mAngleSmear[5][21]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mAngleSmearCos[5][21]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mAngleSmearBeta[5][21]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mAngleSmearCosBeta[5][21]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH3D *h3_mAngleSmearBetaCos[5][21]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt

    TH3D *h3_Random[3];

    TH1D *h_mPhiMassMC[20];
    TH1D *h_mPhiMass[20][42];
    TH2D *h_mPtRes[2][42];

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
    int mStudy;

    int flag_eff;
    int flag_eff_PtEtaPhi;
    int flag_eff_Cos;
    
    int mpt_first;
    int mpt_last;

    double mpt_low[10];
    double mpt_up[10];
 
  ClassDef(StEffHistMangerHelicityGlobal,1)
};

#endif
