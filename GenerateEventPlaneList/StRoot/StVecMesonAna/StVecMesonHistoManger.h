#ifndef StVecMesonHistoManger_h
#define StVecMesonHistoManger_h

#include "StMessMgr.h"
#include <map>

class TH1F;

typedef std::map<TString,TH1F*> TH1FMap;

class StVecMesonHistoManger
{
  public:
    StVecMesonHistoManger(Int_t EP_mode);
    ~StVecMesonHistoManger();

    void InitSys(Int_t X_flag, Int_t mode);
    void FillSys(Float_t pt, Float_t eta, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t Mass2, Double_t reweight, Int_t X_flag, Int_t mode);
    void FillDcaSys(Float_t,Float_t,Int_t);
    void FillSigSys(Float_t,Float_t,Int_t);
    void WriteSys(Int_t X_flag, Int_t mode);

    void InitSys_EP(Int_t X_flag, Int_t mode);
    void FillSys_EP(Float_t pt, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t Mass2, Double_t reweight, Int_t X_flag, Int_t mode);
    void WriteSys_EP(Int_t X_flag, Int_t mode);

  private:
    // spin alignment analysis
    // 0 = pt bin
    // 1 = centrality: 9 = 20%-60%, 0-8 from RefMultCorr 
    // 2 = cos(theta*)
    // 3 = dca => 2.0, 2.5 and 3.0
    // 4 = nSigmaKaon => 2.5, 2.0 and 3.0
    TH1FMap h_mMass2;
    TH1FMap h_mMass2_EP;

    // raw pt spectra
    // 0 = pt bin
    // 1 = pt bin finer
    // 2 = centrality: 9 = 20%-60%, 0-8 from RefMultCorr 
    // 3 = dca => 2.0, 2.5 and 3.0
    // 4 = nSigmaKaon => 2.5, 2.0 and 3.0
    TH1FMap h_mMass_Spec;

    // event plane resolution correction
    // 0 = centrality: 0-8 from RefMultCorr 
    // 1 = dca => 2.0, 2.5 and 3.0
    // 2 = nSigmaKaon => 2.5, 2.0 and 3.0
    TH1FMap h_mMass_Yields;

    TH1FMap h_mDca;
    TH1FMap h_mSigKaon;

    std::string EP_string;

  ClassDef(StVecMesonHistoManger,1)
};
#endif
