#ifndef StVecMesonFillCorr_h
#define StVecMesonFillCorr_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"
#include "TLorentzVector.h"

class TLorentzVector;
class TProfile2D;
class TFile;

class StVecMesonFillCorr : public TObject
{
  public:
    StVecMesonFillCorr(Int_t energy);
    virtual ~StVecMesonFillCorr();

    // ReCenter Correction
    void InitReCenterCorrection();

    TVector2 getReCenterPar_East(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);
    TVector2 getReCenterPar_West(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);
    TVector2 getReCenterPar_Full(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);
 
    TVector2 calq2Vector(TLorentzVector);
    Float_t getWeight(TLorentzVector);

    // Shift Correction
    void InitShiftCorrection();
    bool passTrackEtaNumCut(Int_t, Int_t); // Num of Tracks East, Num of Tracks West with eta_gap = +/- 0.05
    bool passTrackFullNumCut(Int_t, Int_t, Int_t); // Num of Tracks Full, Num of Tracks East, Num of Tracks West

    Float_t AngleShift(Float_t Psi_raw);
    Float_t calShiftAngle2East_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2West_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2Full_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2Full_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign, TLorentzVector lTrack, bool UsedEP);

    // Resolution Correction
    void InitResolutionCorr();
    Float_t getResolution2_EP(Int_t); // centrality
    Float_t getResolution2_Full_EP(Int_t);

  private:
    TFile *mInPutFile_ReCenter; // input file for ReCenter Correction
    TFile *mInPutFile_Shift; // input file for Shift Correction
    TFile *mInPutFile_Res; // input file for Resolution Correction

    Int_t mEnergy;

    static TString mVStr[4];
    static TString mOrder;

  ClassDef(StVecMesonFillCorr,1)
};
#endif
