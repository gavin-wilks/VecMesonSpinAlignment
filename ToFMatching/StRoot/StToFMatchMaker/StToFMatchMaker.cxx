#include "StRoot/StToFMatchMaker/StToFMatchMaker.h"
#include "StRoot/StToFMatchMaker/StToFMatchCut.h"
#include "StRoot/StToFMatchMaker/StToFMatchHistoManger.h"
#include "StRoot/StToFMatchMaker/StUtility.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
//#include "StRoot/StPicoEvent/StPicoV0.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>

ClassImp(StToFMatchMaker)

StRefMultCorr* StToFMatchMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StToFMatchMaker::StToFMatchMaker(const char* name, StPicoDstMaker *picoMaker, const char *jobCounter, const int energy)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mEnergy = energy;

  mOutPut_ToFMatch = Form("file_%s_ToFMatch_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter);
}

//----------------------------------------------------------------------------- 
StToFMatchMaker::~StToFMatchMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StToFMatchMaker::Init() 
{
  mRefMultCorr = new StRefMultCorr("refmult");
 
  mToFMatchCut = new StToFMatchCut(mEnergy);
  mToFMatchHistoManger = new StToFMatchHistoManger();

  mFile_ToFMatch = new TFile(mOutPut_ToFMatch.Data(),"RECREATE");
  mFile_ToFMatch->cd();
  mToFMatchHistoManger->InitQA();
  mToFMatchHistoManger->InitHist();
  mScaleFactor_nSigma = vmsa::mSigScaleMap[mEnergy];

  mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StToFMatchMaker::Finish() 
{
  if(mOutPut_ToFMatch != "")
  {
    mFile_ToFMatch->cd();
    mToFMatchHistoManger->WriteQA();
    mToFMatchHistoManger->WriteHist();
    mFile_ToFMatch->Close();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StToFMatchMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StToFMatchMaker::Make() 
{
  if(!mPicoDstMaker) 
  {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) 
  {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent)
  {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // RefMult
  Int_t runId = mPicoEvent->runId();
  Int_t refMult = mPicoEvent->refMult();
  Float_t vz = mPicoEvent->primaryVertex().z();
  Float_t zdcX = mPicoEvent->ZDCx();
  const unsigned short nBTofMatched = mPicoEvent->nBTOFMatch();

  mRefMultCorr->init(runId);

  if(mRefMultCorr->isBadRun( runId ))
  {
    LOG_ERROR << "Bad Run! Skip!" << endm;
    return kStErr;
  }

  bool isPileUpEvent = false;
  // IMPORTANT: vertex position is needed for Au+Au 19.6 GeV 2019
  if (mRefMultCorr->isPileUpEvent( refMult, nBTofMatched, vz ) ) isPileUpEvent = true;
  mRefMultCorr->initEvent(refMult,vz,zdcX); 

  const Int_t cent9 = mRefMultCorr->getCentralityBin9();       // 0: 70-80%, 1: 60-70%, ..., 6: 10-20%, 7: 5-10%, 8: 0-5%

  // Event Cut
  //const  double refMultCorr = mToFMatchCut->getRefMultReweight(vz, refMult);
  //const Int_t cent9 = mToFMatchCut->getCentrality(refMultCorr);

  // Event Cut
  if(!isPileUpEvent && mToFMatchCut->passEventCut(mPicoDst) && cent9 > -0.5) // event cut
  {

    const int runIndex = mUtility->findRunIndex(runId);
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StUtility! Skip!" << endm;
      return kStErr;
    }

    const Int_t nTracks = mPicoDst->numberOfTracks();
    const Double_t reweight = mRefMultCorr->getWeight();
    //const Int_t nToFMatched = mToFMatchCut->getMatchedToF();

    mToFMatchHistoManger->FillQA_Event(vz,refMult);

    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      float Mass2 = mToFMatchCut->getMass2(track, mPicoDst);
      float dEdx = track->dEdx();
      float p = track->pMom().Mag();
      int polarity = track->charge();
      int charge = 0;
      if(polarity < 0) charge = 1;

      float pt = track->pMom().Perp();
      float eta = track->pMom().PseudoRapidity();
      float phi = track->pMom().Phi();
      if(mToFMatchCut->passTrackBasic(track,mPicoEvent))
      {
	mToFMatchHistoManger->FillQA_Detector(dEdx,Mass2,p*polarity);
      }

      /*
      if(mToFMatchCut->passSigPionCut(track,mScaleFactor_nSigma)) // pion QA
      {
	mToFMatchHistoManger->FillQA_Pion(dEdx,Mass2,p*polarity);
	mToFMatchHistoManger->Fill_TPC(charge,1,cent9,pt,eta,phi);
	if(mToFMatchCut->passToFMatchCut(track))
	{
	  mToFMatchHistoManger->Fill_ToF(charge,1,cent9,pt,eta,phi);
	}
      }
      if(mToFMatchCut->passSigProntonCut(track,mScaleFactor_nSigma)) // proton QA
      {
	mToFMatchHistoManger->FillQA_Proton(dEdx,Mass2,p*polarity);
	mToFMatchHistoManger->Fill_TPC(charge,2,cent9,pt,eta,phi);
	if(mToFMatchCut->passToFMatchCut(track))
	{
	  mToFMatchHistoManger->Fill_ToF(charge,2,cent9,pt,eta,phi);
	}
      }
      */

      if(mToFMatchCut->passSigKaonCut(track,mPicoEvent,mScaleFactor_nSigma)) // Kaon QA
      {
	mToFMatchHistoManger->FillQA_Kaon(dEdx,Mass2,p*polarity);
	mToFMatchHistoManger->Fill_TPC(charge,0,cent9,pt,eta,phi);
	if(mToFMatchCut->passToFMatchCut(track,mPicoEvent,mPicoDst))
	{
	  mToFMatchHistoManger->Fill_ToF(charge,0,cent9,pt,eta,phi);
	}
      }
    }
  }

  return kStOK;
}

