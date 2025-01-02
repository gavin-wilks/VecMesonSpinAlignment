#ifndef StToFMatchHistoManger_h
#define StToFMatchHistoManger_h

#include "TObject.h"
#include "StMessMgr.h"
#include "../Utility/type.h"

class TH1F;
class TH1D;
class TH2F;
class TH3D;
class TH3F;

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

    void InitHistQA();
    void FillEventHistQA(float reweight, int cent, float vx, float vy, float vz, float ntof, float refmult);
    void FillTrackHistQA(float reweight, int cent, int charge, float dca, float nhits, float nhitsratio, float p, float dedx);
    void WriteHistQA();    

    void InitHist();
    int findEta(float eta);
    int findPhi(float phi);
    void Fill_TPC(int charge, int pid, int cent, float pt, float eta, float phi);
    void Fill_ToF(int charge, int pid, int cent, float pt, float eta, float phi);
    void WriteHist();
    
  private:

    // QA plots

    //TH2F *h_mDEdx_Pion[10];
    //TH2F *h_mMass2_Pion[10];
    //TH2F *h_mMass2_Kaon[10];
    //TH2F *h_mDEdx_Proton[10];
    //TH2F *h_mMass2_Proton[10];

    TH1F *h_mDca[2][10];
    TH1F *h_mNHits[2][10];
    TH1F *h_mNHitsRatio[2][10];
    TH2F *h_mDEdx_Kaon[2][10];

    TH3F *h_mVertex[10];
    TH1F *h_mNToFMatch[10];
    TH1F *h_mRefMult[10]; 


    TH3D *h_mTracks_TPC[2][3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality => 9 for miniBias <= temporary
    TH3D *h_mTracks_ToF[2][3][10];
    TH1DMap h_mCounts_TPC; // counts for TPC tracks
    TH1DMap h_mCounts_ToF; // counts for ToF tracks
    TH1D *h_FrameEta_ToF;
    TH1D *h_FramePhi_ToF;

  ClassDef(StToFMatchHistoManger,1)
};
#endif
