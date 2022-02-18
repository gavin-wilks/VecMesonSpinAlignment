#ifndef StVecMesonCut_h
#define StVecMesonCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

class StVecMesonCut : public TObject
{
  public:
    StVecMesonCut(const int energy);
    virtual ~StVecMesonCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*, StPicoEvent*);
    bool passSigPionCut(StPicoTrack*, Float_t);
    bool passSigKaonCut(StPicoTrack*, Float_t);
    bool passSigProntonCut(StPicoTrack*, Float_t);
    bool passTrackMeson(StPicoTrack*, StPicoEvent*);
    bool isMinBias(StPicoEvent*);
    bool isPileUpEvent(int, int);
    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();
    Float_t getBeta(StPicoTrack*, StPicoDst*);
    Float_t getPrimaryMass2(StPicoTrack*, StPicoDst*);
    Float_t getV0Mass2(StPicoTrack*, StPicoDst*);

  private:
    static StRefMultCorr *mRefMultCorr;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;
    Int_t mEnergy;
   
    double pl19[5][5] = {{-22.706,0.863809,-0.00036812,5.97123e-06,-1.27439e-08},   //-145 < vz < -87  cm
                       {-15.5924,0.604207,0.00131806,-2.04754e-06,5.73182e-10},   //-87  < vz < -29  cm 
                       {-13.1158,0.504708,0.00187998,-4.7317e-06,4.82342e-09},    //-29  < vz <  29  cm
                       {-15.9411,0.615061,0.00118242,-1.48902e-06,-2.29371e-10},  // 29  < vz <  87  cm 
                       {-22.468,0.839062,-0.000106213,4.93946e-06,-1.14503e-08}}; // 87  < vz <  145 cm

    double ph19[5][5] = {{33.7733,1.75938,-0.00285868,8.50344e-06,-1.19215e-08},  //-145 < vz < -87  cm
                       {20.5875,1.67371,-0.00307534,7.93756e-06,-8.46257e-09},  //-87  < vz < -29  cm 
                       {15.1015,1.53929,-0.00269203,6.48876e-06,-6.06074e-09},  //-29  < vz <  29  cm
                       {20.7719,1.67316,-0.00315093,8.35824e-06,-9.14822e-09},  // 27  < vz <  87  cm 
                       {33.4926,1.79373,-0.00319461,9.56613e-06,-1.31049e-08}}; // 87  < vz <  145 cm

    //******** POINTWISE VZ CORRECTION *********
    double vz_corr[29];
    vz_corr[0]=0.998171;  // (-145,-135)cm
    vz_corr[1]=0.99364;   // (-135,-125)cm
    vz_corr[2]=0.991578;  // (-125,-115)cm
    vz_corr[3]=0.990252;  // (-115,-105)cm
    vz_corr[4]=0.990494;  // (-105,-95)cm
    vz_corr[5]=0.990065;  // (-95,-85)cm
    vz_corr[6]=0.990332;  // (-85,-75)cm
    vz_corr[7]=0.996478;  // (-75,-65)cm
    vz_corr[8]=0.999687;  // (-65,-55)cm
    vz_corr[9]=0.998645;  // (-55,-45)cm
    vz_corr[10]=0.993835; // (-45,-35)cm
    vz_corr[11]=0.996273; // (-35,-25)cm
    vz_corr[12]=0.998307; // (-25,-15)cm
    vz_corr[13]=0.999295; // (-15,-5)cm
    vz_corr[14]=1;        // (-5,5)cm
    vz_corr[15]=1.00056;  // (5,15)cm
    vz_corr[16]=1.00019;  // (15,25)cm
    vz_corr[17]=0.999894; // (25,35)cm
    vz_corr[18]=0.998907; // (35,45)cm
    vz_corr[19]=1.00468;  // (45,55)cm
    vz_corr[20]=1.0055;   // (55,65)cm
    vz_corr[21]=1.00142;  // (65,75)cm
    vz_corr[22]=0.996247; // (75,85)cm
    vz_corr[23]=0.995789; // (85,95)cm
    vz_corr[24]=0.996513; // (95,105)cm
    vz_corr[25]=0.996928; // (105,115)cm
    vz_corr[26]=0.998196; // (115,125)cm
    vz_corr[27]=1.00097;  // (125,135)cm
    vz_corr[28]=1.00699;  // (135,145)cm

    ClassDef(StVecMesonCut,1)
};
#endif
