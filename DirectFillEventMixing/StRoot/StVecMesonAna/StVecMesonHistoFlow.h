#ifndef StVecMesonHistoFlow_h
#define StVecMesonHistoFlow_h

#include "StMessMgr.h"
#include <map>

class TH1F;

typedef std::map<TString,TH1F*> TH1FMap;

class StVecMesonHistoFlow
{
  public:
    StVecMesonHistoFlow();
    ~StVecMesonHistoFlow();

    void InitSys(Int_t X_flag, Int_t mode);
    void FillSys(Float_t pt, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Int_t nHitSys, Float_t Res2, Float_t Mass2, Double_t reweight, Int_t X_flag, Int_t mode);
    void WriteSys(Int_t X_flag, Int_t mode);

  private:
    // spin alignment analysis
    // 0 = pt bin
    // 1 = centrality: 9 = 20%-60%, 0-8 from RefMultCorr 
    // 2 = cos(theta*)
    // 3 = dca => 2.0, 2.5 and 3.0
    // 4 = nSigmaKaon => 2.5, 2.0 and 3.0
    TH1FMap h_mMass2;

    TH1FMap h_mMass_Yields;

  ClassDef(StVecMesonHistoFlow,1)
};
#endif
