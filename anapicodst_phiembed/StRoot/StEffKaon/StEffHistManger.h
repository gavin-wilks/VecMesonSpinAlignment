#ifndef StEffHistManger_h
#define StEffHistManger_h
#include "TObject.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
//#include "../../../Utility/type.h"
#include <map>
#include <string>

class TH1F;
class TH2F;
class TH3F;
class TH2F;

typedef std::map<std::string,TH1F*> TH1FMap;

class StEffHistManger : public TObject
{
  public:
    StEffHistManger();
    virtual ~StEffHistManger();
    void InitHist();
    float AngleShift(float);
    void FillHistMc(double,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float);
    void FillHistRc(double,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int);
    void FillHistPt(int,float,float,float);
    //void CalEfficiency();
    //void CalEffPtEtaPhi();
    TH1F* CalEffError(TH1F*,TH1F*,std::string);
    void WriteHist();

    void InitHistQA();
    void FillEventHistQA(float reweight, int cent, float vx, float vy, float vz, float ntof, float refmult);
    void FillTrackHistQA(float reweight, int cent, int charge, float dca, float nhits, float nhitsratio, float p, float dedx, int step);
    void WriteHistQA();    

  private:

    TH1F *h_mDca[10][2][10];
    TH1F *h_mNHits[10][2][10];
    TH1F *h_mNHitsRatio[10][2][10];
    TH2F *h_mDEdx_Kaon[10][2][10];

    TH3F *h_mVertex[10];
    TH1F *h_mNToFMatch[10];
    TH1F *h_mRefMult[10]; 


    TH3F *h_mMcTracks[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3F *h_mRcTracks[10][10];
    TH3F *h_mMcTracksKplus[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3F *h_mRcTracksKplus[10][10];
    TH3F *h_mMcTracksKminus[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3F *h_mRcTracksKminus[10][10];

    TH2F *h_mPtGl[10]; // McPt vs. gRcPt
    TH2F *h_mPtPr[10]; // McPt vs. pRcPt

    TH1F *h_mMcEffPt[10]; // pt distritbution as a function of centrality
    TH1F *h_mRcEffPt[10];
    TH1F *h_mEffPt[10];

    TH1F *h_mMcEffEta[10]; // eta distritbution as a function of centrality
    TH1F *h_mRcEffEta[10];
    TH1F *h_mEffEta[10];

    TH1F *h_mMcEffPhi[10]; // phi distritbution as a function of centrality
    TH1F *h_mRcEffPhi[10];
    TH1F *h_mEffPhi[10];

    TH2F *h_mpTv2Weights[10];
    TH3F *h_mpTyv2Weights[10];
    TH2F *h_mpTv2Weights_Var[10][10];

    TH1F *h_FrameEta;
    TH1F *h_FramePhi;
    TH1F *h_FrameEtaKaon;
    TH1F *h_FramePhiKaon;

    TH1FMap h_mMcEffPEP;  
    TH1FMap h_mRcEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi
    TH1FMap h_mEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi

    TH2F *h_mMcEffPtY[10]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
    TH2F *h_mRcEffPtY[10][10];
    TH2F *h_mKaon_MC[2][10];
    TH2F *h_mKaon_RC[10][2][10];

    TH3F *h_mKaon_RC_m2[10][2][10];

    TH2F *h_mMcEffPtPhiPsi[10][10];
    TH3F *h_mMcEffPtYPhiPsi[10][10];

    TH1F *h_mMcEffPhiS[10][6];
    TH1F *h_mRcEffPhiS[10][10][6];
    TH1F *h_mEffPhiS[10][10][6];

    TH1F *h_mMcEffCos[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1F *h_mRcEffCos[10][10][10];
    TH1F *h_mEffCos[10][10][10];

    TH2F *h_mMcEffCosCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH2F *h_mRcEffCosCosH[10][10][10];
    TH2F *h_mEffCosCosH[10][10][10];

    TH1F *h_mMcEffCosH[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1F *h_mRcEffCosH[10][10][10];
    TH1F *h_mEffCosH[10][10][10];

    TH2F *h_mMcEffCosPhiPrime[10][10]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2F *h_mRcEffCosPhiPrime[10][10][10];
    TH2F   *h_mEffCosPhiPrime[10][10][10];

    TH2F *h_mMcEffCosHPhiPrime[10][10]; // helicity efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2F *h_mRcEffCosHPhiPrime[10][10][10];
    TH2F   *h_mEffCosHPhiPrime[10][10][10];

    TH1F *h_mMcEffPhiSY[10][6];
    TH1F *h_mRcEffPhiSY[10][10][6];
    TH1F   *h_mEffPhiSY[10][10][6];

    TH1F *h_mMcEffCosY[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1F *h_mRcEffCosY[10][10][10];
    TH1F   *h_mEffCosY[10][10][10];

    TH2F *h_mMcEffCosCosHY[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH2F *h_mRcEffCosCosHY[10][10][10];
    TH2F   *h_mEffCosCosHY[10][10][10];

    TH1F *h_mMcEffCosHY[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1F *h_mRcEffCosHY[10][10][10];
    TH1F   *h_mEffCosHY[10][10][10];

    TH2F *h_mMcEffCosPhiPrimeY[10][10]; // efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2F *h_mRcEffCosPhiPrimeY[10][10][10];
    TH2F   *h_mEffCosPhiPrimeY[10][10][10];

    TH2F *h_mMcEffCosHPhiPrimeY[10][10]; // helicity efficiency vs CosThetaStar & phiprime as a function of centrality and pt
    TH2F *h_mRcEffCosHPhiPrimeY[10][10][10];
    TH2F   *h_mEffCosHPhiPrimeY[10][10][10];

    int flag_eff;
    int flag_eff_PtEtaPhi;

  ClassDef(StEffHistManger,1)
};

#endif
