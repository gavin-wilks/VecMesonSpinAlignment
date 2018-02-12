#ifndef StZdcSmdProManger_h
#define StZdcSmdProManger_h

#include "TVector2.h"
#include "StMessMgr.h"
#include <string>

class TProfile2D;
class TProfile;

class StZdcSmdProManger
{
  public:
    StZdcSmdProManger();
    virtual ~StZdcSmdProManger();

    // ReCenter Correction
    void InitReCenter();
    void FillReCenterEast(TVector2 qVector, int Cent9, int RunIndex); // vz_sign = vertex pos/neg
    void FillReCenterWest(TVector2 qVector, int Cent9, int RunIndex);
    void WriteReCenter();

    // Shift Correction for East/West
    void InitShift();
    void FillShiftEast(TVector2 qVector, int Cent9, int RunIndex);
    void FillShiftWest(TVector2 qVector, int Cent9, int RunIndex);
    void WriteShift();

    // Shift Correction for Full
    void InitShiftFull();
    void FillShiftFull(TVector2 qVector, int Cent9, int RunIndex);
    void WriteShiftFull();

    // Event Plane Resolution for East/West
    void InitResolution();
    void FillResolution(TVector2 QEast, TVector2 QWest, int Cent9);
    void WriteResolution();

    //  Directed Flow QA
    void InitDirectedFlow();
    void FillDirectedFlow(int Cent9, float eta, float pt, float v1, float resolution, float reweight);
    void WriteDirectedFlow();

  private:
    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mQEastVertical; // 0 = vertex pos/neg
    TProfile2D *p_mQEastHorizontal;
    TProfile2D *p_mQWestVertical;
    TProfile2D *p_mQWestHorizontal;

    // Shift Correction for East/West
    TProfile2D *p_mQEastCos[20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQEastSin[20];
    TProfile2D *p_mQWestCos[20];
    TProfile2D *p_mQWestSin[20];

    // Shift Correction for East/West
    TProfile2D *p_mQFullCos[20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQFullSin[20];

    // event plane resolution for East/West
    TProfile *p_mResolution;

    TProfile *p_mDirectedFlow[9];
    TProfile *p_mDirectedFlowCom;

    static string mVStr[2];

    ClassDef(StZdcSmdProManger,1)
};

#endif
