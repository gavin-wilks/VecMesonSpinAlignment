#include "StRoot/StToFMatchMaker/StToFMatchMaker.h"
#include "StRoot/StToFMatchMaker/StToFMatchCut.h"
#include "StRoot/StToFMatchMaker/StToFMatchHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoV0.h"
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
StToFMatchMaker::StToFMatchMaker(const char* name, StPicoDstMaker *picoMaker, const int jobCounter, const int energy)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mEnergy = energy;

  mOutPut_ToFMatch = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/file_%s_ToFMatch_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
}

//----------------------------------------------------------------------------- 
StToFMatchMaker::~StToFMatchMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StToFMatchMaker::Init() 
{
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mToFMatchCut = new StToFMatchCut(mEnergy);
  mToFMatchHistoManger = new StToFMatchHistoManger();

  mFile_ToFMatch = new TFile(mOutPut_ToFMatch.Data(),"RECREATE");
  mFile_ToFMatch->cd();
  mToFMatchHistoManger->InitQA();
  mToFMatchHistoManger->InitHist();
  mScaleFactor_nSigma = vmsa::mSigScaleMap[mEnergy];

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
  mRefMultCorr->init(runId);
  if(mEnergy == 6) mRefMultCorr->initEvent(refMult,vz,zdcX); // for 200 GeV
  if(mEnergy != 6) mRefMultCorr->initEvent(refMult,vz); // for BES Energy

  // Event Cut
  if(mToFMatchCut->passEventCut(mPicoDst)) // event cut
  {
    const Int_t nTracks = mPicoDst->numberOfTracks();
    const Int_t cent9 = mRefMultCorr->getCentralityBin9();
//    if(cent9 < 0) cout << cent9 << endl;
    const Double_t reweight = mRefMultCorr->getWeight();
    const Int_t nToFMatched = mToFMatchCut->getMatchedToF();

    mToFMatchHistoManger->FillQA_Event(vz,refMult);

    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      float Mass2 = mToFMatchCut->getMass2(track);
      float dEdx = track->dEdx();
      float p = track->pMom().mag();
      int polarity = track->charge();
      int charge = 0;
      if(polarity < 0) charge = 1;

      float pt = track->pMom().perp();
      float eta = track->pMom().pseudoRapidity();
      float phi = track->pMom().phi();
      if(mToFMatchCut->passTrackBasic(track))
      {
	mToFMatchHistoManger->FillQA_Detector(dEdx,Mass2,p*polarity);
      }

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

      if(mToFMatchCut->passSigKaonCut(track,mScaleFactor_nSigma)) // Kaon QA
      {
	mToFMatchHistoManger->FillQA_Kaon(dEdx,Mass2,p*polarity);
	mToFMatchHistoManger->Fill_TPC(charge,0,cent9,pt,eta,phi);
	if(mToFMatchCut->passToFMatchCut(track))
	{
	  mToFMatchHistoManger->Fill_ToF(charge,0,cent9,pt,eta,phi);
	}
      }
    }
  }

  return kStOK;
}

