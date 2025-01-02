#ifndef StEpdEpManager_h
#define StEpdEpManager_h

#include "TObject.h"
#include "TVector2.h"
#include "StEventPlaneCons.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
class TNtuple;
class TProfile;

class StEpdEpManager : public TObject
{
  public:
    StEpdEpManager(int energy);
    virtual ~StEpdEpManager();
    void clearEpdEp();
    void initEpdEp(int Cent9, int RunIndex);

    void readResolution();
    double getResSub(int Cent9, int order);
    double getResSubErr(int Cent9, int order);

  private:
    //Event Plane method
    int mEnergy;
    int mCent9;
    int mRunIndex;

    // EP resolution
    TFile *mInPutFile_Res;
    double mEpdSubResVal[9][vmsa::mEpdEpOrder];
    double mEpdSubResErr[9][vmsa::mEpdEpOrder];

    TProfile *p_mEpdSubRes[vmsa::mEpdEpOrder];
 
  ClassDef(StEpdEpManager,1)
};

#endif
