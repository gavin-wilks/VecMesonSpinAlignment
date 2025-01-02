#ifndef StVecMesonHistoManger_h
#define StVecMesonHistoManger_h

#include "StMessMgr.h"
#include <map>
#include <string>

class TH1F;
class TH2F;


typedef std::map<std::string,TH2F*> TH2FMap;
typedef std::map<std::string,TH1F*> TH1FMap;


class StVecMesonHistoManger
{
  public:
    StVecMesonHistoManger();
    virtual ~StVecMesonHistoManger();

    void InitQA();
    void InitPID();
    void FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p);
    void FillQA_Event(Float_t vz, Float_t refMult);
    void FillPID(int cent, short charge, double pt, double eta, double phi, double nsig, double mass2);
    void FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West);
    void FillEP_Full(Float_t Psi2_Full);
    void WriteQA();
    void WritePID();

    void InitEP();
    void FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift);
    void FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift);
    void WriteEP();
    
  private:


    // QA plots
    TH2F *h_mDEdx;
    TH2F *h_mMass2;

    //TH2F *h_mPID[2][20][10][12]; // first is the charge 0 --> -1 , 1 --> 1 | Then the pT bin | Then the pesudorapidity bin
    TH2FMap h_mPID; // first is the charge 0 --> -1 , 1 --> 1 | Then the pT bin | Then the pesudorapidity bin
    TH1FMap h_mPID_1D; // first is the charge 0 --> -1 , 1 --> 1 | Then the pT bin | Then the pesudorapidity bin

    TH1F *h_mEastRaw;
    TH1F *h_mWestRaw;
    TH1F *h_mFullRaw;
    TH1F *h_mVz;
    TH1F *h_mRefMult;

    // event plane distribution
    TH1F *h_mEastReCenter;
    TH1F *h_mWestReCenter;
    TH1F *h_mRanAReCenter;
    TH1F *h_mRanBReCenter;
    TH1F *h_mFullReCenter;

    TH1F *h_mEastShift;
    TH1F *h_mWestShift;
    TH1F *h_mRanAShift;
    TH1F *h_mRanBShift;
    TH1F *h_mFullShift;

  ClassDef(StVecMesonHistoManger,1)
};
#endif
