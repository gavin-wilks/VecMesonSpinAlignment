#ifndef StEffHistManger_h
#define StEffHistManger_h
#include "StMessMgr.h"
#include "TObject.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/type.h"
#include <map>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"

class TProfile;
class TH1D;
class TH2D;
class TH3D;
class TObject;

class StEffHistManger : public TObject
{
  public:
    StEffHistManger(int energy, int pid, int mode, int startpt, int stoppt);
    virtual ~StEffHistManger();
    void InitHist();
    void FillHistMc(int,float,float,float,float,float,float,float);
    void FillHistRc(int,float,float,float,float,float,float,float);
    //float AngleShift(float);
    void CalEffCosThetaStar();
    TH1D* CalEffError(TH1D*,TH1D*,std::string);
    void WriteHist();

  private:
    TH1D *h_mMcEffCosSPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TH1D *h_mRcEffCosSPt[10][10];
    TH1D   *h_mEffCosSPt[10][10];
    TProfile *p_mMcCos2BetaPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mMcCos2BetaPPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mRcCos2BetaPt[10][10];
    TProfile *p_mRcCos2BetaPPt[10][10];
    TProfile *p_mMcCos4BetaPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mMcCos4BetaPPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mRcCos4BetaPt[10][10];
    TProfile *p_mRcCos4BetaPPt[10][10];
    TProfile *p_mMcCos2BetaCos4BetaPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mMcCos2BetaPCos4BetaPPt[10][10]; // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mRcCos2BetaCos4BetaPt[10][10];
    TProfile *p_mRcCos2BetaPCos4BetaPPt[10][10];
    TH1D *h_mMcEffCosSPt_RP[10][10]; // efficiency vs CosThetaStar w.r.t. RP as a function of centrality and pt
    TH1D *h_mRcEffCosSPt_RP[10][10];
    TH1D   *h_mEffCosSPt_RP[10][10];
    TGraphAsymmErrors *g_mEffCosSPt[10][10];

    TH1D *h_mMcEffCosPt[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosPt[10][10];
    TH1D   *h_mEffCosPt[10][10];
    TGraphAsymmErrors *g_mEffCosPt[10][10];

    TH1D *h_mMcEffCosSY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TH1D *h_mRcEffCosSY[10][10][vmsa::y_total];
    TH1D   *h_mEffCosSY[10][10][vmsa::y_total];
    TProfile *p_mMcCos2BetaY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mMcCos2BetaPY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mRcCos2BetaY[10][10][vmsa::y_total];
    TProfile *p_mRcCos2BetaPY[10][10][vmsa::y_total];
    TProfile *p_mMcCos4BetaY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mMcCos4BetaPY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mRcCos4BetaY[10][10][vmsa::y_total];
    TProfile *p_mRcCos4BetaPY[10][10][vmsa::y_total];
    TProfile *p_mMcCos2BetaCos4BetaY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mMcCos2BetaPCos4BetaPY[10][10][vmsa::y_total];    // efficiency vs CosThetaStar w.r.t. EP as a function of centrality and pt
    TProfile *p_mRcCos2BetaCos4BetaY[10][10][vmsa::y_total];
    TProfile *p_mRcCos2BetaPCos4BetaPY[10][10][vmsa::y_total];
    TH1D *h_mMcEffCosSY_RP[10][10][vmsa::y_total]; // efficiency vs CosThetaStar w.r.t. RP as a function of centrality and pt
    TH1D *h_mRcEffCosSY_RP[10][10][vmsa::y_total];
    TH1D   *h_mEffCosSY_RP[10][10][vmsa::y_total];
    TGraphAsymmErrors *g_mEffCosSY[10][10][vmsa::y_total];

    TH1D *h_mMcEffCosY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosY[10][10][vmsa::y_total];
    TH1D   *h_mEffCosY[10][10][vmsa::y_total];
    TGraphAsymmErrors *g_mEffCosY[10][10][vmsa::y_total];

    int mEnergy, mMode;
    int flag_eff_Cos;

    int mStartPt;
    int mStopPt;    

    int mpt_first;
    int mpt_last;

    float mpt_low[10];
    float mpt_up[10];
 
  ClassDef(StEffHistManger,1)
};

#endif
