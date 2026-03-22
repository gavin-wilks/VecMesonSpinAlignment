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
    void FillHistMc(int,float,float,float,float);
    void FillHistRc(int,float,float,float,float);
    //float AngleShift(float);
    void CalEffCosThetaStar();
    TH1D* CalEffError(TH1D*,TH1D*,std::string);
    void WriteHist();

  private:
    TH1D *h_mMcEffCosSPt[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosSPt[10][10];
    TH1D   *h_mEffCosSPt[10][10];

    TH1D *h_mMcEffCosPt[10][10]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosPt[10][10];
    TH1D   *h_mEffCosPt[10][10];

    TH1D *h_mMcEffCosSY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosSY[10][10][vmsa::y_total];
    TH1D   *h_mEffCosSY[10][10][vmsa::y_total];

    TH1D *h_mMcEffCosY[10][10][vmsa::y_total]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCosY[10][10][vmsa::y_total];
    TH1D   *h_mEffCosY[10][10][vmsa::y_total];

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
