#ifndef StToFMatchHistoManger_h
#define StToFMatchHistoManger_h

#include "TObject.h"
#include "StMessMgr.h"
#include "../Utility/type.h"

class TH1F;
class TH1D;
class TH2F;
class TH3D;

class StToFMatchHistoManger : public TObject
{
  public:
    StToFMatchHistoManger();
    virtual ~StToFMatchHistoManger();

    void InitQA();
    void FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p);
    void FillQA_Pion(Float_t dEdx, Float_t Mass2, Float_t p);
    void FillQA_Kaon(Float_t dEdx, Float_t Mass2, Float_t p);
    void FillQA_Proton(Float_t dEdx, Float_t Mass2, Float_t p);
    void FillQA_Event(Float_t vz, Float_t refMult);
    void WriteQA();
    
    void InitHist();
    int findEta(float eta);
    int findPhi(float phi);
    void Fill_TPC(int charge, int pid, int cent, float pt, float eta, float phi);
    void Fill_ToF(int charge, int pid, int cent, float pt, float eta, float phi);
    void WriteHist();
    
  private:

    // QA plots
    TH2F *h_mDEdx;
    TH2F *h_mMass2;

    TH2F *h_mDEdx_Pion;
    TH2F *h_mMass2_Pion;
    TH2F *h_mDEdx_Kaon;
    TH2F *h_mMass2_Kaon;
    TH2F *h_mDEdx_Proton;
    TH2F *h_mMass2_Proton;

    TH1F *h_mVz;
    TH1F *h_mRefMult;

    TH3D *h_mTracks_TPC[2][3][9]; // pt, eta, phi distribution as a function of charge | pid | centrality
    TH3D *h_mTracks_ToF[2][3][9];
    TH1DMap h_mCounts_TPC; // counts for TPC tracks
    TH1DMap h_mCounts_ToF; // counts for ToF tracks
    TH1D *h_FrameEta;
    TH1D *h_FramePhi;

  ClassDef(StToFMatchHistoManger,1)
};
#endif
