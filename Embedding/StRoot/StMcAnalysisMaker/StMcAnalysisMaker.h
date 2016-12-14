#ifndef ST_MCANALYSISMAKER_H
#define ST_MCANALYSISMAKER_H

#include <vector>
#include "TString.h"

#include "StChain/StMaker.h"

class TFile;
class TNtuple;
class TH3F;

class StMcTrack;
class StTrack;
class StGlobalTrack;
class StAssociationMaker;
class StMcEvent;
class StEvent;
class StMuDst;
class StDedxPidTraits;
class StRefMultCorr;
class StPrimaryTrack;

class StMcAnalysisMaker : public StMaker
{
private:
   TString mOutfileName;
   StRefMultCorr* mRefMultCorrUtil;
   StMuDst*       mMuDst;
   std::vector<float> firedTriggersIndices;
   double mField; //.. magnetic field
   int    mCentrality;
   bool mFillTpcHitsNtuple;

   TFile* mFile;
   TNtuple* mTracks;
   TNtuple* mEventCount; //.. For counting purposes
   TNtuple* mTpcNtuple;

   TH3F* hTpcHitsDiffXVsPadrowVsSector;                                                                                                                                                                  
   TH3F* hTpcHitsDiffYVsPadrowVsSector;
   TH3F* hTpcHitsDiffZVsPadrowVsSector;

   StMcEvent* mMcEvent;
   StEvent* mEvent;
   StAssociationMaker* mAssoc;

   StTrack const* findPartner(StMcTrack*, int&) const;
   StMcTrack const* findPartner(StGlobalTrack*, int&) const;
   StDedxPidTraits const* findDedxPidTraits(StTrack const*) const;
   int getNHitsDedx(StTrack const*) const;

   bool passTrigger();
   int  fillEventCounts(float nRTracks = -1, float nMcTracks = -1);
   int  fillTracks(int& nRTracks, int& nMcTracks);
   void fillMcTrack(float* array,int& idx,StMcTrack const*);
   void fillRcTrack(float* array,int& idx,StTrack const* gTrk,StPrimaryTrack const* pTrk,int const ncom);
   unsigned int  getHftTruth(StMcTrack const*,StTrack const*) const;
   void getDca(StTrack const*,float& dca, float& dcaXY, float& dcaZ) const;

   bool isGoodMcTrack(StMcTrack const*) const;

public:
   StMcAnalysisMaker (const char *name="StMcAnalysisMaker", const char *title="event/StMcAnalysisMaker");

   int Init();
   int Make();
   int Finish();

   void setOutFileName(std::string);
   void fillTpcHitsNtuple(bool t=false);
   void setRefMultCorr(StRefMultCorr*);

   ClassDef(StMcAnalysisMaker, 0)
};

inline void StMcAnalysisMaker::setOutFileName(std::string s){ mOutfileName = s.c_str();}
inline void StMcAnalysisMaker::setRefMultCorr(StRefMultCorr* refmultCorrUtil){ mRefMultCorrUtil = refmultCorrUtil;}
#endif
