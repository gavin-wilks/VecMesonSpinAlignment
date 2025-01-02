#ifndef StVecMesonHistoManger_h
#define StVecMesonHistoManger_h

#include "StMessMgr.h"
#include "../Utility/phi_data_constants_19GeV.h"
#include <map>

class TH1F;
class TH2F;
class TH3F;
class TProfile;

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TH2F*> TH2FMap;
typedef std::map<TString,TProfile*> TProMap;

class StVecMesonHistoManger
{
  public:
    StVecMesonHistoManger(Int_t EP_mode);
    ~StVecMesonHistoManger();

    void InitSys(Int_t X_flag, Int_t mode);
    void InitPhiSys(Int_t X_flag, Int_t mode);
    void FillSysKplus(Int_t cent9, Float_t phi_pt, Float_t Cos2PhiStarPhi, Float_t kaon_pt, Float_t kaon_phi, Float_t kaon_y, Float_t kaon_eta, Float_t InvMass, Double_t reweight, Int_t X_flag);
    void FillSysKminus(Int_t cent9, Float_t phi_pt, Float_t Cos2PhiStarPhi, Float_t kaon_pt, Float_t kaon_phi, Float_t kaon_y, Float_t kaon_eta, Float_t InvMass, Double_t reweight, Int_t X_flag);
    //void FillSys(Int_t cent9, Float_t phi_pt, Float_t Cos2PhiStarPhi, Float_t kaon_pt, Float_t kaon_phi, Float_t kaon_y, Float_t kaon_eta, Float_t InvMass, Double_t reweight, Int_t X_flag);
    void FillSys(Int_t cent9, Float_t phi_pt, Float_t phi_y, Float_t Cos2PhiStarPhi, Float_t kaon_pt, Float_t kaon_phi, Float_t kaon_y, Float_t kaon_eta, Int_t kcharge, Float_t InvMass, Double_t reweight, Int_t X_flag);
    void FillPhiSys(Int_t cent9, Float_t pt, Float_t phi, Float_t y, Float_t phi_star, Float_t Cos2PhiStarPhi, Float_t InvMass, Double_t reweight, Int_t X_flag);
    void FillDeltaSys(Int_t cent9, Float_t phi_pt, Float_t Cos2PhiStarPhi, Float_t phi, Float_t eta, Float_t kaonp_pt, Float_t kaonp_phi, Float_t kaonp_eta, Float_t kaonm_pt, Float_t kaonm_phi, Float_t kaonm_eta, Float_t InvMass, Double_t reweight, Int_t X_flag);
    void FillDcaSys(Float_t,Float_t,Int_t);
    void FillSigSys(Float_t,Float_t,Int_t);
    void WritePhiSys(Int_t X_flag, Int_t mode);
    void WriteSys(Int_t X_flag, Int_t mode);

    void InitSys_EP(Int_t X_flag, Int_t mode);
    void FillSys_EP(Float_t pt, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t Mass2, Double_t reweight, Int_t X_flag, Int_t mode);
    void WriteSys_EP(Int_t X_flag, Int_t mode);

    void InitHistQA(Int_t, Int_t);
    void FillEventHistQA(float reweight, int cent, float vx, float vy, float vz, float ntof, float refmult);
    void FillTrackHistQA(float phi_pt, float phi_y, float InvMass, float reweight, int cent, int charge, float dca, float nhits, float nhitsratio, float p, float dedx, int X_flag);
    void WriteHistQA();    

  private:
    // spin alignment analysis
    // 0 = pt bin
    // 1 = centrality: 9 = 20%-60%, 0-8 from RefMultCorr 
    // 2 = cos(theta*)
    // 3 = dca => 2.0, 2.5 and 3.0
    // 4 = nSigmaKaon => 2.5, 2.0 and 3.0
    TH1F *h_mDcaTrack[2][10];
    TH1F *h_mNHits[2][10];
    TH1F *h_mNHitsRatio[2][10];
    TH2F *h_mDEdx_Kaon[2][10];

    TH3F *h_mVertex[10];
    TH1F *h_mNToFMatch[10];
    TH1F *h_mRefMult[10]; 

    TProMap p_mMass2;
    TH1FMap h_mMass2;
    TH1FMap h_mMass2_EP;

    //TH2F *h_mKaon[4][10][data_constants::kaon_pt_bins];
    TH2FMap h_mKaon;
    
    TH2FMap h_mPhi;
    TH1FMap h_mPhi1D;

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

    std::string EP_string = "";

  ClassDef(StVecMesonHistoManger,1)
};
#endif
