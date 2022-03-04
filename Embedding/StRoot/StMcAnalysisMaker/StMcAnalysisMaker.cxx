#include "TFile.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TSystem.h"

#include "StarClassLibrary/StParticleDefinition.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StParticleTypes.hh"
#include "StBFChain/StBFChain.h"

#include "StEvent/StEventTypes.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StTrackTopologyMap.h"
#include "StEvent/StEventSummary.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StTpcDedxPidAlgorithm.h"
#include "StEventUtilities/StuRefMult.hh"

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"

#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcContainers.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcEvent.hh"

#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"

#include "StRefMultCorr/StRefMultCorr.h"
#include "StMcAnaCuts.h"
#include "StMcAnalysisMaker.h"

ClassImp(StMcAnalysisMaker);

StMcAnalysisMaker::StMcAnalysisMaker(const char *name, const char *title): StMaker(name),
   mRefMultCorrUtil(NULL), mPicoDst(NULL), mField(-999), mCentrality(-999), 
   mFile(NULL), mTracks(NULL), mEventCount(NULL), mMcEvent(NULL), mEvent(NULL), mAssoc(NULL)
{
   LOG_INFO << "StMcAnalysisMaker() - DONE" << endm;
}
//__________________________________
int StMcAnalysisMaker::Init()
{
   if (!mOutfileName.Length())
   {
      // StBFChain* bfChain = (StBFChain *) StMaker::GetChain();
      //
      // if (!bfChain) return kStFatal;
      //
      // mOutfileName = bfChain->GetFileIn();

      if (mOutfileName.Length())
      {
         LOG_INFO << mOutfileName << endm;
         mOutfileName = gSystem->BaseName(mOutfileName.Data());
         mOutfileName = mOutfileName.ReplaceAll(".event.root", "");
         mOutfileName = mOutfileName.ReplaceAll(".geant.root", "");
         mOutfileName = mOutfileName.ReplaceAll(".PicoDst.root", "");
      }
      else
      {
         mOutfileName = "mcAnalysis";
      }
   }

   mOutfileName = mOutfileName.ReplaceAll(".root", "");
   mFile = new TFile(Form("%s.McAna.root", mOutfileName.Data()), "recreate");
   assert(mFile && !mFile->IsZombie());

   mAssoc = (StAssociationMaker*)GetMaker("StAssociationMaker");
   if (!mAssoc)
   {
      LOG_ERROR << "Could not get StAssociationMaker" << endm;
      return kStErr;
   }

   for (int ii = 0; ii < McAnaCuts::maxNumberOfTriggers; ++ii)
   {
      firedTriggersIndices.push_back(-1);
   };

   const char* evtVarList = "runId:eventId:mcVx:mcVy:mcVz:vx:vy:vz:vzVpd:"
     "centrality:gRefMult:RefMult:posRefMult:negRefMult:zdc:bbc:nMcTracks:nRcTracks:magField:t0:t1:t2:t3:t4:t5";

   const char* varlist = "McPt:McP:McEta:McY:McPhi:McGeantId:eventGenLabel:StartVtxX:StartVtxY:StartVtxZ:StopVtxX:StopVtxY:StopVtxZ:" // MC Kaon
     "gRcPt:gRcEta:gRcPhi:gRcNfit:gRcNmax:gRcNcom:gRcNdedx:gRcDedx:gRcNsigKP:gRcNsigKM:gRcDca:gRcDcaXY:gRcDcaZ:" // global RC K+
     "pRcPt:pRcEta:pRcPhi:pRcNfit:pRcNmax:pRcNcom:pRcNdedx:pRcDedx:pRcNsigKP:pRcNsigKM:pRcDca:pRcDcaXY:pRcDcaZ"; // primary RC K+

   mEventCount = new TNtuple("eventCount","eventCount",evtVarList);

   mTracks = new TNtuple("Kaon","",varlist);

   LOG_INFO << "Init() - DONE" << endm;

   return kStOk;
}

//__________________________________
int StMcAnalysisMaker::Make()
{
   StPicoDstMaker* picoDstMaker = (StPicoDstMaker*)GetMaker("PicoDst");

   if (!picoDstMaker)
   {
      LOG_WARN << " No PicoDstMaker, will try to take all event information from StEvent" << endm;
      mPicoDst = NULL;
   }
   else
   {
     mPicoDst = (StPicoDst*)picoDstMaker->picoDst();
   }

   if(!mPicoDst || !mPicoDst->event())
   {
     LOG_WARN << "PicoDst or mPicoDst->event() is missing, will try to take all event information from StEvent" << endm;
     mPicoDst = NULL;
   }

   mMcEvent = (StMcEvent*)GetDataSet("StMcEvent");

   if (!mMcEvent)
   {
      LOG_WARN << "No StMcEvent" << endm;
      return kStWarn;
   }

   mEvent = (StEvent*)GetDataSet("StEvent");
   if (!mEvent)
   {
      LOG_WARN << "No StEvent" << endm;
      return kStWarn;
   }

   mField = (float)mEvent->summary()->magneticField();

   if(mRefMultCorrUtil && mPicoDst)
   {
     mRefMultCorrUtil->init(mEvent->runId());

     mRefMultCorrUtil->initEvent(mPicoDst->event()->refMult(),
     // mRefMultCorrUtil->initEvent(mPicoDst->event()->grefmult(),
                                  mEvent->primaryVertex()->position().z(), 
                                  mEvent->runInfo()->zdcCoincidenceRate());

     mCentrality  = mRefMultCorrUtil->getCentralityBin9();

     if (mRefMultCorrUtil->isBadRun(mEvent->runId()))
     {
       LOG_INFO << "This is a bad run from mRefMultCorrUtil! Skip! " << endm;

       return kStSkip;
     }
   }
   else
   {
     mCentrality = -999;
   }

   // Fill
   int nRcTracks = -1;
   int nMcTracks = -1;

   int fillTracksStatus = kStOk;
   // fillTracksStatus = fillTracks(nRcTracks, nMcTracks);
   if (passTrigger())
   {
     fillTracksStatus = fillTracks(nRcTracks, nMcTracks);
   }
   else
   {
      LOG_INFO << "No interesting triggers. Counting event then skipping." << endm;
   }

   int fillEventCountStatus = fillEventCounts((float)nRcTracks, (float)nMcTracks);

   return fillTracksStatus && fillEventCountStatus;
}

int StMcAnalysisMaker::fillTracks(int& nRcTracks, int& nMcTracks)
{
   nRcTracks = 0;
   nMcTracks = 0;

   LOG_INFO << "Filling " << mMcEvent->tracks().size() << " mcTracks..." << "\n";

   for (unsigned int iTrk = 0;  iTrk < mMcEvent->tracks().size(); ++iTrk)
   {
      StMcTrack* const mcTrack = mMcEvent->tracks()[iTrk];

      if (!mcTrack)
      {
         LOG_WARN << "Empty mcTrack container" << endm;
         continue;
      }

      if (!isGoodMcTrack(mcTrack)) continue;
      ++nMcTracks;

      int nCommonHits = -999;
      StTrack const* const gTrk = findPartner(mcTrack, nCommonHits);
      StPrimaryTrack* rcTrack = 0;
      if(gTrk) rcTrack = (StPrimaryTrack*)gTrk->node()->track(primary); // primary


      float array[220];
      for (int ii = 0; ii < 220; ++ii) array[ii] = -999;

      int idx = 0;

      fillMcTrack(array, idx, mcTrack);

      if (gTrk)
      {
         ++nRcTracks;
         fillRcTrack(array, idx, gTrk, rcTrack, nCommonHits);
      }

      mTracks->Fill(array);
   }

   LOG_INFO << endm;

   return kStOk;
}

void StMcAnalysisMaker::fillMcTrack(float* array, int& idx, StMcTrack const* const mcTrk)
{
  array[idx++] = mcTrk->pt();
  array[idx++] = mcTrk->momentum().mag();
  array[idx++] = mcTrk->pseudoRapidity();
  array[idx++] = mcTrk->rapidity();
  array[idx++] = mcTrk->momentum().phi();
  array[idx++] = mcTrk->geantId();
  array[idx++] = mcTrk->eventGenLabel();
  array[idx++] = mcTrk->startVertex()->position().x();
  array[idx++] = mcTrk->startVertex()->position().y();
  array[idx++] = mcTrk->startVertex()->position().z();
  if (mcTrk->stopVertex())
  {
    array[idx++] = mcTrk->stopVertex()->position().x();
    array[idx++] = mcTrk->stopVertex()->position().y();
    array[idx++] = mcTrk->stopVertex()->position().z();
  }
  else
  {
    idx += 3;
  }
}

void StMcAnalysisMaker::fillRcTrack(float* array, int& idx, StTrack const* const gTrk, StPrimaryTrack const* const rcTrack, int const ncom)
{
  // global RcTracks
  array[idx++] = gTrk ? gTrk->geometry()->momentum().perp() : -999.;
  array[idx++] = gTrk ? gTrk->geometry()->momentum().pseudoRapidity() : -999.;
  array[idx++] = gTrk ? gTrk->geometry()->momentum().phi() : -999.;
  array[idx++] = gTrk ? gTrk->fitTraits().numberOfFitPoints(kTpcId) : -999.;
  array[idx++] = gTrk ? gTrk->numberOfPossiblePoints(kTpcId) : -999.;
  array[idx++] = gTrk ? ncom : -999.;

  // dedx info
  float nDedxPts = -999.;
  float dedx = -999.;
  float nSigKP = -999.;
  float nSigKM = -999.;
  static StTpcDedxPidAlgorithm aplus(McAnaCuts::dedxMethod);
  static StKaonPlus* KPlus = StKaonPlus::instance();
  static StKaonMinus* KMinus = StKaonMinus::instance();
  StParticleDefinition const* prtcl = 0;
  if(gTrk) prtcl = gTrk->pidTraits(aplus);
  if (prtcl)
  {
    nDedxPts = aplus.traits()->numberOfPoints();
    dedx = aplus.traits()->mean();
    nSigKP = aplus.numberOfSigma(KPlus);
    nSigKM = aplus.numberOfSigma(KMinus);
  }

  // array[idx++] = getNHitsDedx(gTrk);
  array[idx++] = nDedxPts;
  array[idx++] = dedx;
  array[idx++] = nSigKP;
  array[idx++] = nSigKM;

  float dca = -999.;
  float dcaXY = -999.;
  float dcaZ = -999.;

  if(gTrk) getDca(gTrk, dca, dcaXY, dcaZ);

  array[idx++] = dca;
  array[idx++] = dcaXY;
  array[idx++] = dcaZ;

  // primary RcTracks
  array[idx++] = rcTrack ? rcTrack->geometry()->momentum().perp() : -999.;
  array[idx++] = rcTrack ? rcTrack->geometry()->momentum().pseudoRapidity() : -999.;
  array[idx++] = rcTrack ? rcTrack->geometry()->momentum().phi() : -999.;
  array[idx++] = rcTrack ? rcTrack->fitTraits().numberOfFitPoints(kTpcId) : -999.;
  array[idx++] = rcTrack ? rcTrack->numberOfPossiblePoints(kTpcId) : -999.;
  array[idx++] = rcTrack ? ncom : -999.;

  // dedx info
  float nDedxPtsPr = -999.;
  float dedxPr = -999.;
  float nSigKPPr = -999.;
  float nSigKMPr = -999.;
  static StTpcDedxPidAlgorithm aplusPr(McAnaCuts::dedxMethod);
  static StKaonPlus* KPlusPr = StKaonPlus::instance();
  static StKaonMinus* KMinusPr = StKaonMinus::instance();
  StParticleDefinition const* prtclPr = 0;
  if(rcTrack) prtclPr = rcTrack->pidTraits(aplusPr);
  if (prtclPr)
  {
    nDedxPtsPr = aplusPr.traits()->numberOfPoints();
    dedxPr = aplusPr.traits()->mean();
    nSigKPPr = aplusPr.numberOfSigma(KPlusPr);
    nSigKMPr = aplusPr.numberOfSigma(KMinusPr);
  }

  // array[idx++] = getNHitsDedx(rcTrack);
  array[idx++] = nDedxPtsPr;
  array[idx++] = dedxPr;
  array[idx++] = nSigKPPr;
  array[idx++] = nSigKMPr;

  float dcaPr = -999.;
  float dcaXYPr = -999.;
  float dcaZPr = -999.;

  if(rcTrack) getDca(rcTrack, dcaPr, dcaXYPr, dcaZPr);

  array[idx++] = dcaPr;
  array[idx++] = dcaXYPr;
  array[idx++] = dcaZPr;
}

bool StMcAnalysisMaker::isGoodMcTrack(StMcTrack const* const mcTrack) const
{
   return mcTrack->geantId() == McAnaCuts::geantId && mcTrack->startVertex()->position().perp() < McAnaCuts::mcTrackStartVtxR;
}

int StMcAnalysisMaker::fillEventCounts(float nRcTracks, float nMcTracks)
{
   float vars[50];

   float vpdVz = -999;
   StBTofHeader* tofheader = 0;
   if (mEvent->btofCollection())  tofheader = mEvent->btofCollection()->tofHeader();
   if (tofheader) vpdVz = tofheader->vpdVz();

   int iVar = 0;
   vars[iVar++] = (float)mEvent->runId();
   vars[iVar++] = (float)mEvent->id();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().z();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().z();
   vars[iVar++] = vpdVz;
   vars[iVar++] = mCentrality;
   vars[iVar++] = mPicoDst? mPicoDst->event()->grefmult() : -999;
   vars[iVar++] = mPicoDst? mPicoDst->event()->refMult() : -999;
   vars[iVar++] = (float)uncorrectedNumberOfPositivePrimaries(*mEvent);
   vars[iVar++] = (float)uncorrectedNumberOfNegativePrimaries(*mEvent);
   vars[iVar++] = (float)mEvent->runInfo()->zdcCoincidenceRate();
   vars[iVar++] = (float)mEvent->runInfo()->bbcCoincidenceRate();
   vars[iVar++] = nMcTracks;
   vars[iVar++] = nRcTracks;
   vars[iVar++] = (float)mEvent->summary()->magneticField() / 10;
   vars[iVar++] = firedTriggersIndices.at(0);
   vars[iVar++] = firedTriggersIndices.at(1);
   vars[iVar++] = firedTriggersIndices.at(2);
   vars[iVar++] = firedTriggersIndices.at(3);
   vars[iVar++] = firedTriggersIndices.at(4);
   vars[iVar++] = firedTriggersIndices.at(5);

   mEventCount->Fill(vars);

   return kStOk;
}

bool StMcAnalysisMaker::passTrigger()
{
   LOG_INFO << "Checking triggers..." << endm;
   bool interesting_event = false;

   if (!mEvent)
   {
      LOG_FATAL << "mEvent doesn't exist" << endm;
   }

   if (McAnaCuts::interesting_triggers.size() == 0)
   {
      LOG_WARN << "No triggers in McAnaCuts::interesting_triggers ... accepting event anyway" << endm;
      return true;
   }

   const StTriggerId* st_trgid = mEvent->triggerIdCollection()->nominal();

   for (unsigned int ii = 0; ii < firedTriggersIndices.size(); ++ii)
   {
      firedTriggersIndices[ii] = -1;
   }

   // Fill interesting triggers
   LOG_INFO << "Interesting fired triggers: " << "\n";

   for (unsigned int ii = 0; ii < st_trgid->maxTriggerIds(); ++ii)
   {
      unsigned int id = st_trgid->triggerId(ii);

      int trgIndex = -1;

      for (unsigned int jj = 0; jj < McAnaCuts::interesting_triggers.size(); ++jj)
      {
         if (McAnaCuts::interesting_triggers[jj] == id)
         {
            trgIndex = jj;
            interesting_event = true;
            LOG_INFO << id << " ";
            break;
         }
      }

      if (trgIndex != -1) firedTriggersIndices[trgIndex] = 1.0;
   }

   LOG_INFO << endm;

   return interesting_event;
}

StTrack const* StMcAnalysisMaker::findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits) const
{
   //..StMcTrack find partner from the StTracks
   pair<mcTrackMapIter, mcTrackMapIter> p = mAssoc->mcTrackMap()->equal_range(mcTrack);

   const StTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (mcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();
      const StTrack* track = k->second->partnerTrack()->node()->track(global);//should be global
      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

void StMcAnalysisMaker::getDca(StTrack const* const rcTrack, float& dca, float& dcaXY, float& dcaZ) const
{
   StPhysicalHelixD helix = rcTrack->geometry()->helix();
   dca = helix.distance(mEvent->primaryVertex()->position());
   dcaXY = helix.geometricSignedDistance(mEvent->primaryVertex()->position().x(), mEvent->primaryVertex()->position().y());

   StThreeVectorF dcaPoint = helix.at(helix.pathLength(mEvent->primaryVertex()->position()));
   dcaZ = dcaPoint.z() - mEvent->primaryVertex()->position().z();
}

unsigned int StMcAnalysisMaker::getHftTruth(StMcTrack const* const mcTrack, StTrack const* const rcTrack) const
{
   bool istTruth = true;
   bool pxlTruth1 = true;
   bool pxlTruth2 = true;

   StPtrVecHit rcIstHits = rcTrack->detectorInfo()->hits(kIstId);
   int nRcIstHits = (int)rcIstHits.size();
   for (int iHit = 0; iHit < nRcIstHits; ++iHit)
   {
      if (rcIstHits[iHit]->idTruth() != mcTrack->key())
      {
         istTruth = false;
         break;
      }
   }

   StPtrVecHit rcPxlHits = rcTrack->detectorInfo()->hits(kPxlId);
   int nRcPxlHits = (int)rcPxlHits.size();
   for (int iHit = 0; iHit < nRcPxlHits; ++iHit)
   {
      if (rcPxlHits[iHit]->idTruth() != mcTrack->key())
      {
         StThreeVectorF pos = rcPxlHits[iHit]->position();

         float const R = pow(pos.x(), 2.0) + pow(pos.y(), 2.0);
         if (R > 3.5 * 3.5) pxlTruth2 = false;
         else pxlTruth1 = false;
      }
   }

   unsigned int hftTruth = 0;
   if (pxlTruth1) hftTruth |= (1 << 0);
   if (pxlTruth2) hftTruth |= (1 << 1);
   if (istTruth)  hftTruth |= (1 << 2);

   return hftTruth;
}

StMcTrack const* StMcAnalysisMaker::findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits) const
{
   //.. StGlobalTracks find partner from StMcTracks.
   //.. See example from StRoot/StMcAnalysisMaker
   pair<rcTrackMapIter, rcTrackMapIter> p = mAssoc->rcTrackMap()->equal_range(rcTrack);

   const StMcTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (rcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();

      const StMcTrack* track = k->second->partnerMcTrack();

      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

int StMcAnalysisMaker::Finish()
{
   mFile->cd();
   mFile->Write();
   mFile->Close();
   return kStOk;
}

StDedxPidTraits const* StMcAnalysisMaker::findDedxPidTraits(StTrack const* const rcTrack) const
{
   StDedxPidTraits* pid = 0;
   StPtrVecTrackPidTraits traits = rcTrack->pidTraits(kTpcId);

   for (unsigned int ii = 0; ii < traits.size(); ++ii)
   {
      pid = dynamic_cast<StDedxPidTraits*>(traits[ii]);
      if (pid && pid->method() == McAnaCuts::dedxMethod) break;
   }

   return pid;
}

int StMcAnalysisMaker::getNHitsDedx(StTrack const* const t) const
{
   int ndedx = -999;
   StPtrVecTrackPidTraits pidTraits = t->pidTraits(kTpcId);

   if (pidTraits.size())
   {
      StDedxPidTraits* pid;
      for (unsigned int ii = 0; ii < pidTraits.size(); ++ii)
      {
         pid = dynamic_cast<StDedxPidTraits*>(pidTraits[ii]);

         if (pid && (pid->method() == McAnaCuts::dedxMethod))
         {
            ndedx = pid->numberOfPoints();            //number of dedx hits
            break;
         }
      }
   }

   return ndedx;
}
