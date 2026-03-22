#include "StRoot/StVecMesonMaker/StVecMesonMaker.h"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/StVecMesonMaker/StVecMesonProManger.h"
#include "StRoot/StVecMesonMaker/StVecMesonCorr.h"
#include "StRoot/StVecMesonMaker/StVecMesonHistoManger.h"
#include "StRoot/StVecMesonMaker/StVecMesonTree.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
//#include "StRoot/StPicoEvent/StPicoV0.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
//#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StRoot/StVecMesonMaker/StUtility.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>

ClassImp(StVecMesonMaker)

StRefMultCorr* StVecMesonMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StVecMesonMaker::StVecMesonMaker(const char* name, StPicoDstMaker *picoMaker, const char *jobCounter, const Int_t Mode, const Int_t energy, const Int_t flag_ME, const Int_t flag_PID)
  : StMaker(name)
{
  cout << "StVecMesonMaker::StVecMesonMaker" << endl;
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;
  mFlag_PID = flag_PID;  

  if(mMode == 0)
  {
    mOutPut_ReCenterPar = Form("file_%s_ReCenterPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 1)
  {
    mOutPut_ShiftPar = Form("file_%s_ShiftPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
  }
  if(mMode == 2)
  {
    mOutPut_Resolution = Form("file_%s_Resolution_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
  }
  if(mMode == 3)
  { 
    mOutPut_PID = Form("file_%s_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[mFlag_PID].c_str(),vmsa::MixEvent[mFlag_ME].Data(),jobCounter); 
  }
  if(mMode == 4)
  { 
    mOutPut_PID = Form("file_PID_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[mFlag_PID].c_str(),jobCounter); 
  }
}

//----------------------------------------------------------------------------- 
StVecMesonMaker::~StVecMesonMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StVecMesonMaker::Init() 
{
  cout << "Start Initialization" << endl;
  mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index 
  cout << "mUtility->initRunIndex()" << endl;
  mRefMultCorr = new StRefMultCorr("refmult"); 
  //if(!mRefMultCorr)
  //{
  //  mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  //}
  cout << "StRefMultCorr" << endl;
  mVecMesonCut = new StVecMesonCut(mEnergy);
  cout << "StVecMesonCut" << endl;
  mVecMesonProManger = new StVecMesonProManger();
  cout << "StVecMesonProManger" << endl;
  mVecMesonCorrection = new StVecMesonCorrection(mEnergy);
  cout << "StVecMesonCorrection" << endl;
  mVecMesonHistoManger = new StVecMesonHistoManger();
  cout << "StVecMesonHistoManger" << endl;

  if(mMode == 0)
  {
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.Data(),"RECREATE");
    mFile_ReCenterPar->cd();
    mVecMesonProManger->InitReCenter();
    mVecMesonHistoManger->InitQA();
  }

  if(mMode == 1)
  {
    //mUsedTrackCounter = 0;
    mVecMesonCorrection->InitReCenterCorrection();
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.Data(),"RECREATE");
    mVecMesonProManger->InitShift();
  }

  if(mMode == 2)
  {
    mUsedTrackCounter = 0;
    mVecMesonCorrection->InitReCenterCorrection();
    mVecMesonCorrection->InitShiftCorrection();
    mVecMesonProManger->InitResolution();
    mVecMesonHistoManger->InitEP();
    mFile_Resolution = new TFile(mOutPut_Resolution.Data(),"RECREATE");
  }

  if(mMode == 3)
  {
    mVecMesonTree = new StVecMesonTree(mEnergy);
    cout << "StVecMesonTree" << endl;
    mFile_PID = new TFile(mOutPut_PID.Data(),"RECREATE");
    mFile_PID->cd();
    if(mFlag_PID == 0)
    {
      mVecMesonTree->InitPhi();
      cout << "InitPhi" << endl;
    }
    mVecMesonCorrection->InitReCenterCorrection();
    cout << "InitReCenterCorrection" << endl;
    mVecMesonCorrection->InitShiftCorrection();
    cout << "InitShiftCorrection" << endl;
    mVecMesonCorrection->InitResolutionCorr();
    cout << "InitResolutionCorr" << endl;
  }
  if(mMode == 4)
  {
    mFile_PID = new TFile(mOutPut_PID.Data(),"RECREATE");
    mFile_PID->cd();
    mVecMesonCorrection->InitReCenterCorrection();
    cout << "InitReCenterCorrection" << endl;
    mVecMesonCorrection->InitShiftCorrection();
    cout << "InitShiftCorrection" << endl;
    mVecMesonCorrection->InitResolutionCorr();
    cout << "InitResolutionCorr" << endl;
    mVecMesonHistoManger->InitPID();
    cout << "InitPID" << endl;
  }

  cout << "Finished Initialization" << endl;

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StVecMesonMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();
      mVecMesonProManger->WriteReCenter();
      mVecMesonHistoManger->WriteQA();
      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_ShiftPar != "")
    {
      mFile_ShiftPar->cd();
      mVecMesonProManger->WriteShift();
      mFile_ShiftPar->Close();
    }
  }
  if(mMode == 2)
  {
    if(mOutPut_Resolution != "")
    {
      mFile_Resolution->cd();
      mVecMesonHistoManger->WriteEP();
      mVecMesonProManger->WriteResolution();
      mFile_Resolution->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_PID != "")
    {
      mFile_PID->cd();
      if(mFlag_PID == 0) mVecMesonTree->WriteMass2Phi();
      mFile_PID->Close();
    }
  }
  if(mMode == 4)
  {
    if(mOutPut_PID != "")
    {
      mFile_PID->cd();
      mVecMesonHistoManger->WritePID();
      mFile_PID->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StVecMesonMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StVecMesonMaker::Make() 
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
  Float_t vx = mPicoEvent->primaryVertex().x();
  Float_t vy = mPicoEvent->primaryVertex().y();
  Float_t vz = mPicoEvent->primaryVertex().z();
  Float_t zdcX = mPicoEvent->ZDCx();
  const unsigned short nBTofMatched = mPicoEvent->nBTOFMatch();
  mRefMultCorr->init(runId);
  //mRefMultCorr->initEvent(refMult,vz,zdcX); 

  // vz sign
  //Int_t vz_sign;
  //if(vz > 0.0)
  //{
  //  vz_sign = 0;
  //}
  //else
  //{
  //  vz_sign = 1;
  //}

  // runIndex
  //mRunIdEventsDb = StRunIdEventsDb::Instance((Float_t)mPicoEvent->energy(),(Float_t)mPicoEvent->year());
  //const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(runId); // expensive
//  cout << runIndex << endl;
//  cout << mRunIdEventsDb->getTotalNrRunIds() << endl;

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
  //const  double refMultCorr = mVecMesonCut->getRefMultReweight(vz, refMult);
  //const Int_t cent9 = mVecMesonCut->getCentrality(refMultCorr);
  if(!isPileUpEvent && mVecMesonCut->passEventCut(mPicoDst) && cent9 > -0.5) // event cut
  {
    int vz_sign = 0;
    if(vz > -70.0 && vz <= -30.0) vz_sign = 0;
    if(vz > -30.0 && vz <= 0.0  ) vz_sign = 1;
    if(vz > 0.0   && vz <= +30.0) vz_sign = 2;
    if(vz < +70.0 && vz >  +30.0) vz_sign = 3;

    const Double_t reweight = mRefMultCorr->getWeight();           // Retrieve weight

    //cout << "passed event cut" << endl;
    const int runIndex = mUtility->findRunIndex(runId); // find run index for a specific run
    
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StUtility! Skip!" << endm;
      return kStErr;
    }

    const Int_t nTracks = mPicoDst->numberOfTracks();
    //    if(cent9 < 0) cout << cent9 << endl;
    //const Double_t reweight = mVecMesonCut->getEventWeight(cent9, refMultCorr);
    //const Int_t nToFMatched = mVecMesonCut->getMatchedToF();
    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      if(mMode == 0)
      {
	Float_t eta = track->pMom().PseudoRapidity();
	if(fabs(eta) < vmsa::mEtaMax && track->gDCA(vx,vy,vz) < 3.0)
	{
	  Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track, mPicoDst);
	  Float_t dEdx = track->dEdx();
	  Float_t p = track->pMom().Mag();
	  mVecMesonHistoManger->FillQA_Detector(dEdx,Mass2,p);
	}
      }
      if(mVecMesonCut->passTrackEP(track,mPicoEvent)) // track cut
      {
	if(mMode == 0) // fill re-center parameter
	{
	  Float_t pt = track->pMom().Perp();

	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    TVector2 q2Vector_Full = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManger->FillTrackFull(q2Vector_Full,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_FullRaw(track,cent9,runIndex);
	  }

	  if(mVecMesonCorrection->passTrackEtaEast(track)) // neg eta sub
	  {
	    TVector2 q2Vector_East = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManger->FillTrackEast(q2Vector_East,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_EastRaw(track,cent9,runIndex);
	  }
	  if(mVecMesonCorrection->passTrackEtaWest(track)) // pos eta sub
	  {
	    TVector2 q2Vector_West = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManger->FillTrackWest(q2Vector_West,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_WestRaw(track,cent9,runIndex);
	  }
	}
	else // calculate Q Vector after recentering for full event and eta sub event
	{
	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    mVecMesonCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
	    mUsedTrackCounter++;
	  }
	  if(mVecMesonCorrection->passTrackEtaEast(track)) // neg eta sub
	  {
	    mVecMesonCorrection->addTrack_East(track,cent9,runIndex,vz_sign);
	  }
	  if(mVecMesonCorrection->passTrackEtaWest(track)) // pos eta sub
	  {
	    mVecMesonCorrection->addTrack_West(track,cent9,runIndex,vz_sign);
	  }
	}
      }
    }

    if(mMode == 0) // fill raw EP
    {
      mVecMesonHistoManger->FillQA_Event(vz,refMult);
      if(mVecMesonCorrection->passTrackEtaNumRawCut())
      {
	TVector2 Q2East = mVecMesonCorrection->getQVectorRaw(0); // 0 = eta_gap, 1 = east/west
	Float_t Psi2_East = 0.5*TMath::ATan2(Q2East.Y(),Q2East.X());
	TVector2 Q2West = mVecMesonCorrection->getQVectorRaw(1); // 0 = eta_gap, 1 = east/west
	Float_t Psi2_West = 0.5*TMath::ATan2(Q2West.Y(),Q2West.X());
	mVecMesonHistoManger->FillEP_Eta(Psi2_East,Psi2_West);
      }
      if(mVecMesonCorrection->passTrackFullNumRawCut())
      {
	TVector2 Q2Full = mVecMesonCorrection->getQVectorRaw(2);
	Float_t Psi2_Full = 0.5*TMath::ATan2(Q2Full.Y(),Q2Full.X());
	mVecMesonHistoManger->FillEP_Full(Psi2_Full);
      }
    }

    if(mMode == 1)
    {
      // full event shift parameter
      if(mVecMesonCorrection->passTrackFullNumCut())
      {
	for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
	{
	  TVector2 Psi2Vector_Full_EP = mVecMesonCorrection->calPsi2_Full_EP(k);
	  mVecMesonProManger->FillEventFull_EP(Psi2Vector_Full_EP,cent9,runIndex,vz_sign,k);
	}
      }

      // eta sub event shift parameter
      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	for(Int_t k = 0; k < 5; k++)
	{
	  TVector2 Psi2Vector_East_EP = mVecMesonCorrection->calPsi2_East_EP(k);
	  mVecMesonProManger->FillEventEast_EP(Psi2Vector_East_EP,cent9,runIndex,vz_sign,k);

	  TVector2 Psi2Vector_West_EP = mVecMesonCorrection->calPsi2_West_EP(k);
	  mVecMesonProManger->FillEventWest_EP(Psi2Vector_West_EP,cent9,runIndex,vz_sign,k);
	}
      }
    }

    if(mMode == 2) // calculate resolution for eta_sub and random sub event plane
    {
      // calculate Q vector after recentering for Random Sub Event
      //cout << "Got to the make section" << endl;
      //cout << mUsedTrackCounter << endl;
      Int_t iTrack[mUsedTrackCounter];
      Float_t ranCounter = (Float_t)mUsedTrackCounter/2.0 - 1;
      //cout << "Loop over mUsedTrackCounter" << endl;
      for(Int_t i = 0; i < mUsedTrackCounter; i++)
      {
        iTrack[i] = i;
      }
      //cout << "Randomly shuffle" << endl;
      std::srand(time(0));
      std::random_shuffle(iTrack,iTrack+mUsedTrackCounter);

      mUsedTrackCounter = 0;
      //cout << "Loop over tracks" << endl;
      for(Int_t i = 0; i < nTracks; i++) // track loop
      {
	StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	if(mVecMesonCut->passTrackEP(track,mPicoEvent)) // track cut
	{
          //cout << "passtrackEP" << endl;
	  if(mVecMesonCorrection->passTrackFull(track))
	  {
            //cout << "pass track full" << endl;
	    if((Float_t)iTrack[mUsedTrackCounter] > ranCounter) // Sub Event A
	    {
              //cout << "Sub event A" << endl;
	      mVecMesonCorrection->addTrack_A(track,cent9,runIndex,vz_sign);
	    }
	    else // Sub Event B
	    {
              //cout << "Sub event B" << endl;
	      mVecMesonCorrection->addTrack_B(track,cent9,runIndex,vz_sign);
	    }
	    mUsedTrackCounter++;
	  }
	}
      }
      mUsedTrackCounter = 0;
      //cout << "Calculte resolution" << endl;
      // calculate resolution
      TVector2 QVecEast = mVecMesonCorrection->getQVector(0);
      Float_t Psi2East_ReCenter = 0.5*TMath::ATan2(QVecEast.Y(),QVecEast.X());
      Float_t Psi2East_Shift = mVecMesonCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign);

      TVector2 QVecWest = mVecMesonCorrection->getQVector(1);
      Float_t Psi2West_ReCenter = 0.5*TMath::ATan2(QVecWest.Y(),QVecWest.X());
      Float_t Psi2West_Shift = mVecMesonCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign);

      TVector2 QVecFull = mVecMesonCorrection->getQVector(2);
      Float_t Psi2Full_ReCenter = 0.5*TMath::ATan2(QVecFull.Y(),QVecFull.X());
      Float_t Psi2Full_Shift = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);

      TVector2 QVecRanA = mVecMesonCorrection->getQVector(3);
      Float_t Psi2RanA_ReCenter = 0.5*TMath::ATan2(QVecRanA.Y(),QVecRanA.X());
      Float_t Psi2RanA_Shift = mVecMesonCorrection->calShiftAngle2A_EP(runIndex,cent9,vz_sign);

      TVector2 QVecRanB = mVecMesonCorrection->getQVector(4);
      Float_t Psi2RanB_ReCenter = 0.5*TMath::ATan2(QVecRanB.Y(),QVecRanB.X());
      Float_t Psi2RanB_Shift = mVecMesonCorrection->calShiftAngle2B_EP(runIndex,cent9,vz_sign);

      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
        //cout << "pass eta num cut" << endl;
	mVecMesonHistoManger->FillEP_Sub(Psi2East_ReCenter,Psi2East_Shift,Psi2West_ReCenter,Psi2West_Shift);
	mVecMesonProManger->FillRes_Sub(cent9,Psi2East_Shift,Psi2West_Shift);
      }

      if(mVecMesonCorrection->passTrackFullNumCut())
      {
        //cout << "pass track full num cut" << endl;
	mVecMesonHistoManger->FillEP_Ran(Psi2RanA_ReCenter,Psi2RanA_Shift,Psi2RanB_ReCenter,Psi2RanB_Shift,Psi2Full_ReCenter,Psi2Full_Shift);
	mVecMesonProManger->FillRes_Ran(cent9,Psi2RanA_Shift,Psi2RanB_Shift);
      }
    }

    if(mMode == 3)
    { // phi meson
      if(mVecMesonCorrection->passTrackFullNumCut() && !(cent9 >= 2 && cent9 <= 5))
      {
	// get QVector of sub event
	TVector2 Q2East = mVecMesonCorrection->getQVector(0); // east
	TVector2 Q2West = mVecMesonCorrection->getQVector(1); // west
	TVector2 Q2Full = mVecMesonCorrection->getQVector(2); // full 
	Int_t NumTrackEast = mVecMesonCorrection->getNumTrack(0);
	Int_t NumTrackWest = mVecMesonCorrection->getNumTrack(1);
	Int_t NumTrackFull = mVecMesonCorrection->getNumTrack(2);
	Int_t NumTrackFullEast = mVecMesonCorrection->getNumTrack(3);
	Int_t NumTrackFullWest = mVecMesonCorrection->getNumTrack(4);

	Float_t Psi2 = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);
	// get N_prim, N_non_prim, N_Tof_match
	Int_t N_prim = nBTofMatched;
	Int_t N_non_prim = mPicoEvent->numberOfPrimaryTracks();
	Int_t N_Tof_match = mPicoEvent->numberOfGlobalTracks();

	// pass the event information to StVecMesonTree
	mVecMesonTree->clearEvent();
	mVecMesonTree->passEvent(N_prim, N_non_prim, N_Tof_match);
	// pass re-centered event plane to StVecMesonTree
	mVecMesonTree->passEventPlane(Q2East,Q2West,Q2Full);

	// pass NumOfTrack to StVecMesonTree
	mVecMesonTree->passNumTrack(NumTrackEast,NumTrackWest,NumTrackFull,NumTrackFullEast,NumTrackFullWest);
        if(mFlag_PID == 0)
        {
	  mVecMesonTree->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi2);
        }
      }
    }
    if(mMode == 4)
    { // phi meson
      if(mVecMesonCorrection->passTrackFullNumCut())
      {
        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
          StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
 
          if(mVecMesonCut->passTrackMeson(track,mPicoEvent,0) && track->gDCA(vx,vy,vz) < 2.0)
          {         
            double pT = track->pMom().Perp();
            double phi = track->pMom().Phi();
            if(phi > TMath::Pi())  phi -= (TMath::Pi()*2.0);
            if(phi < -TMath::Pi()) phi += (TMath::Pi()*2.0);
            double eta  = track->pMom().PseudoRapidity();
            double nsig = track->nSigmaKaon();
            double mass2 = mVecMesonCut->getPrimaryMass2(track, mPicoDst);
            short charge = track->charge();

            mVecMesonHistoManger->FillPID(cent9,charge,pT,eta,phi,nsig,mass2);
          }
        }  
      }
    }
    mVecMesonCorrection->clear();
  }
  return kStOK;
}

