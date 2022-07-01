#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdMaker.h"
#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdCut.h"
#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdProManger.h"
#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdCorr.h"
#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdHistoManger.h"
#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdTree.h"
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

#include "StPicoEvent/StPicoEpdHit.h"

#include "StEpdUtil/StEpdEpFinder.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StEpdUtil/StEpdGeom.h"

#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneProManager.h"
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"

#include "StPicoEvent/StPicoEpdHit.h"

ClassImp(StVecMesonEpdMaker)

StRefMultCorr* StVecMesonEpdMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StVecMesonEpdMaker::StVecMesonEpdMaker(const char* name, StPicoDstMaker *picoMaker, const char *jobCounter, const Int_t Mode, const Int_t energy, const Int_t flag_ME, const Int_t flag_PID)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;
  mFlag_PID = flag_PID;  

  if(mMode == 0)
  {
    mInPut_EpdCorrections = "";
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
  }
  if(mMode == 1)
  {
    mInPut_EpdCorrections = Form("StRoot/Utility/EpdCorrections/file_%s_EpdCorrections_0.root",vmsa::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
  }  
  if(mMode == 2)
  {
    mInPut_EpdCorrections = Form("StRoot/Utility/EpdCorrections/file_%s_EpdCorrections_1.root",vmsa::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
    mOutPut_EpdResults = Form("file_%s_EpdResEta_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);  
  }
  if(mMode == 3)
  {
    mInPut_EpdCorrections = Form("StRoot/Utility/EpdCorrections/file_%s_EpdCorrections_0.root",vmsa::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
    mOutPut_EpdResults = Form("file_%s_EpdFlowEta_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);  
  }
  if(mMode == 4)
  {
    mInPut_EpdCorrections = Form("StRoot/Utility/EpdCorrections/file_%s_EpdCorrections_3.root",vmsa::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
    mOutPut_EpdResults = Form("file_%s_EpdFlow_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);  
  }
  if(mMode == 5)
  {
    mInPut_EpdCorrections = Form("StRoot/Utility/EpdCorrections/file_%s_EpdCorrections_0.root",vmsa::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
    mOutPut_EpdResults = Form("file_%s_EpdFlowEta_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);  
  }
  if(mMode == 6)
  {
    mInPut_EpdCorrections = Form("StRoot/Utility/EpdCorrections/file_%s_EpdCorrections_5.root",vmsa::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);
    mOutPut_EpdResults = Form("file_%s_EpdFlow_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mMode,jobCounter);  
  }
 
  //if(mMode == 0)
  //{
  //  mOutPut_ReCenterPar = Form("file_%s_ReCenterPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  //}
  //if(mMode == 1)
  //{
  //  mOutPut_ShiftPar = Form("file_%s_ShiftPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
  //}
  //if(mMode == 2)
  //{
  //  mOutPut_Resolution = Form("file_%s_Resolution_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
  //}
  //if(mMode == 3)
  //{ 
  //  mOutPut_PID = Form("file_%s_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[mFlag_PID].c_str(),vmsa::MixEvent[mFlag_ME].Data(),jobCounter); 
  //}
}

//----------------------------------------------------------------------------- 
StVecMesonEpdMaker::~StVecMesonEpdMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StVecMesonEpdMaker::Init() 
{
  mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index 
  mRefMultCorr = new StRefMultCorr("refmult"); 
  mVecMesonCut = new StVecMesonEpdCut(mEnergy);
  mVecMesonProManger = new StVecMesonEpdProManger();
  mVecMesonCorrection = new StVecMesonEpdCorrection(mEnergy);
  mVecMesonHistoManger = new StVecMesonEpdHistoManger();
  mEventPlaneHistoManager = new StEventPlaneHistoManager();
  mEventPlaneProManager = new StEventPlaneProManager();
  if (mMode < 7 && mMode >= 0) mEpdEpManager = new StEpdEpManager(mEnergy); // initialize EPD EP Manager

  if (mMode < 7 && mMode >= 0)
  {
    mEpFinder = new StEpdEpFinder(10,mOutPut_EpdCorrections.c_str(),mInPut_EpdCorrections.c_str()); // Set the nEventType to 1, since we do not have separation for centrality yet 
    mEpFinder->SetnMipThreshold(0.3);  // recommended by EPD group
    mEpFinder->SetMaxTileWeight(2.0);    // recommended by EPD group
    mEpFinder->SetEpdHitFormat(2);           // 2=pico

    mEpdGeom = new StEpdGeom;
  }

  // On second pass
  //double mResCorr[9] = {0.224327, 0.291191, 0.365294, 0.424694, 0.456142, 0.453143, 0.391865, 0.280205, 0.157329};
  //double mEtaLin[9] = {-0.00455916, -0.00580228, -0.00706776, -0.0075328, -0.00726617, -0.00634251, -0.00438439, -0.00214997, -0.000634274};
  //double mEtaCub[9] = {0.000418717, 0.000568683, 0.000739452, 0.000859493, 0.000899388, 0.000842957, 0.000625815, 0.000327728, 0.000105771};

  double mResCorr[9] = {0.170333, 0.236634, 0.311019, 0.373951, 0.412054, 0.415775, 0.363042, 0.260374, 0.145185};
  double mEtaLin[9] = {-0.00331161, -0.00466038, -0.00593718, -0.0066493, -0.00662886, -0.00585393, -0.00412081, -0.00205545, -0.000625279};
  double mEtaCub[9] = {0.00030976, 0.000457207, 0.000623666, 0.000755677, 0.000815535, 0.000775468, 0.000583732, 0.000308166, 9.95871e-05};  

  for(unsigned int i = 0; i < 9; i++){
    mEtaLin[i] = mEtaLin[i]/mResCorr[i];
    mEtaCub[i] = mEtaCub[i]/mResCorr[i];
  }
  
  if(mMode == 0 || mMode == 1 || mMode == 2)
  {
    double FirstRingWeight[16] = {1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double* PtrRingWeight = FirstRingWeight;
    mEpFinder->SetRingWeights(1, PtrRingWeight);
  }
 
  if(mMode == 2)
  {
    mFile_EpdResults = new TFile(mOutPut_EpdResults.c_str(),"RECREATE");
 
    mEventPlaneProManager->initEpdRes();
    mEventPlaneProManager->initEpdFlowEta();
    mEventPlaneHistoManager->initEpdEp();
    mEventPlaneHistoManager->initEpdQ();
  }

  if(mMode == 3 || mMode == 5)
  {
    TH2D wt1("Order1etaWeight","Order1etaWeight",100,1.5,6.5,10,0,10);
    for (int ix=1; ix<101; ix++){
      for (int iy=1; iy<11; iy++){
        double eta = wt1.GetXaxis()->GetBinCenter(ix);
        if(iy  < 10) wt1.SetBinContent(ix,iy,mEtaLin[iy-1]*eta+mEtaCub[iy-1]*eta*eta*eta);
        if(iy == 10) wt1.SetBinContent(ix,iy,1.0);
      }
    }
    mEpFinder->SetEtaWeights(1,wt1);
  }

  if(mMode == 4 || mMode == 6)
  {
    mFile_EpdResults = new TFile(mOutPut_EpdResults.c_str(),"RECREATE");

    mEventPlaneProManager->initEpdRes();
    mEventPlaneProManager->initEpdFlow();
    mEventPlaneProManager->initEpdFlowEta();
    mEventPlaneHistoManager->initEpdEp();
    mEventPlaneHistoManager->initEpdQ();
    mVecMesonProManger->InitD12();

    TH2D wt1("Order1etaWeight","Order1etaWeight",100,1.5,6.5,10,0,10);
    for (int ix=1; ix<101; ix++){
      for (int iy=1; iy<11; iy++){
        double eta = wt1.GetXaxis()->GetBinCenter(ix);
        if(iy  < 10) wt1.SetBinContent(ix,iy,mEtaLin[iy-1]*eta+mEtaCub[iy-1]*eta*eta*eta);
        if(iy == 10) wt1.SetBinContent(ix,iy,1.0);
      }
    }
    mEpFinder->SetEtaWeights(1,wt1);


    mVecMesonTree = new StVecMesonEpdTree(mEnergy);
    //mFile_PID = new TFile(mOutPut_PID.Data(),"RECREATE");
    //mFile_PID->cd();
    if(mFlag_PID == 0)
    {
      mVecMesonTree->InitPhi();
    }
    if(mFlag_PID == 1)
    {
      mVecMesonTree->InitRho();
    }
    if(mFlag_PID == 2)
    {
      mVecMesonTree->InitKStar();
    }
    mVecMesonCorrection->InitReCenterCorrection();
    mVecMesonCorrection->InitShiftCorrection();
    mVecMesonCorrection->InitResolutionCorr();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StVecMesonEpdMaker::Finish() 
{
  if (mMode < 7 && mMode >= 0) mEpFinder->Finish();

  if(mMode == 2)
  {
    mFile_EpdResults->cd();
    mEventPlaneProManager->writeEpdRes();
    mEventPlaneProManager->writeEpdFlowEta();
    mFile_EpdResults->Close();
  }
  if(mMode == 4 || mMode == 6)
  {
    mFile_EpdResults->cd();
    //mEventPlaneProManager->writeEpdRes();
    //mEventPlaneProManager->writeEpdFlow();
    //mEventPlaneProManager->writeEpdFlowEta();
    mEventPlaneHistoManager->writeEpdEp();
    mEventPlaneHistoManager->writeEpdQ();
    mVecMesonProManger->WriteD12();
    //if(mFlag_PID == 0) mVecMesonTree->WriteMass2Phi();
    //if(mFlag_PID == 1) mVecMesonTree->WriteMass2Rho();
    //if(mFlag_PID == 2) mVecMesonTree->WriteMass2KStar();

    mFile_EpdResults->Close();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StVecMesonEpdMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StVecMesonEpdMaker::Make() 
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
  const TVector3 pv = mPicoEvent->primaryVertex();
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

    StEpdEpInfo result;
    if (mMode < 7 && mMode >= 0)
    {
      result = mEpFinder->Results(mPicoDst->picoArray(StPicoArrays::EpdHit),pv,cent9);
      mEpdEpManager->initEpdEp(cent9,runIndex);
    }
    //mMode == 0 sets phi-weighting and mMode == 1 sets the psi shift corrections
    if(mMode == 2 || mMode == 4 || mMode == 6)
    {// Fill TProfiles for EP resolution and numerator for flow calculations.
     // These two TProfile plots will be used to set eta weights in an external macro.
      //cout << "We are in mMode == 2" << endl;

    //  mEventPlaneProManager->fillEpdRes(result,cent9,runIndex);

      if(mPicoDst->numberOfEpdHits() > 0)
      {
        //cout << "Number of Epd Hits = " << mPicoDst->numberOfEpdHits() << endl;
        /*for(unsigned int iEpd = 0; iEpd < mPicoDst->numberOfEpdHits(); iEpd++)
        {

          //cout << "Before grabbing the epdHit" << endl;

          StPicoEpdHit *EpdHit = mPicoDst->epdHit(iEpd);
          if (! EpdHit) continue;

          //cout << "We have an EPD hit" << endl;          

          TVector3 pos = mEpdGeom->RandomPointOnTile(EpdHit->id());
          TVector3 StraightLine = pos - pv;
          double phi = StraightLine.Phi();
          while(phi < 0.0) phi += 2.0*TMath::Pi();
          double eta = StraightLine.Eta();
          if (!(fabs(eta) > 0) || (fabs(eta) > 1000)) continue;

          //cout << "We have found the position of the hit" << endl;

          int ew = (EpdHit->id() < 0)? 0 : 1;  //is EPD east or west

          double nMip = EpdHit->nMIP(); // Weight of the track
          double nMipEff = nMip; // Put a max on nMip in calculation to prevent slow tracks from overcontributing 
          if(nMipEff > 3.0) nMipEff = 3.0;
          else if(nMipEff < 0.3) continue;

          //cout << "Determined the weight of the track" << endl;

          double v1, v2, v3;
          if (1 == ew){  //west
            v1 = TMath::Cos(1.0*(phi-result.EastPhiWeightedAndShiftedPsi(1)));
            //v2 = TMath::Cos(2.0*(phi-result.EastPhiWeightedAndShiftedPsi(2)));
            //v3 = TMath::Cos(3.0*(phi-result.EastPhiWeightedAndShiftedPsi(3)));

            mEventPlaneProManager->fillEpdFlowEta(eta,v1,cent9,1,nMipEff);
            //mEventPlaneProManager->fillEpdFlowEta(eta,v2,cent9,2,nMipEff); 
            //mEventPlaneProManager->fillEpdFlowEta(eta,v3,cent9,3,nMipEff);

            if(mMode == 4 || mMode == 6)mEventPlaneProManager->fillEpdFlow(v1,cent9,1,runIndex);
            //cout << "Filling west flow" << endl;
          }
          else if(0 == ew){ //east
            v1 = TMath::Cos(1.0*(phi-result.WestPhiWeightedAndShiftedPsi(1)));
            //v2 = TMath::Cos(2.0*(phi-result.WestPhiWeightedAndShiftedPsi(2)));
            //v3 = TMath::Cos(3.0*(phi-result.WestPhiWeightedAndShiftedPsi(3)));

            mEventPlaneProManager->fillEpdFlowEta(-eta,-v1,cent9,1,nMipEff);
            //mEventPlaneProManager->fillEpdFlowEta(eta,v2,cent9,2,nMipEff); 
            //mEventPlaneProManager->fillEpdFlowEta(eta,v3,cent9,3,nMipEff);
            //cout << "Filling east flow" << endl;
            if(mMode == 4 || mMode == 6)mEventPlaneProManager->fillEpdFlow(v1,cent9,1,runIndex);
          }
        }*/
        //cout << "Outside of the loop" << endl;
        mEventPlaneHistoManager->fillEpdEp(result, cent9);
        //cout << "Filled EPD EP" << endl; 
        mEventPlaneHistoManager->fillEpdQ(result, cent9);
        //cout << "Filled EPD Q" << endl;
        //cout << "Successfully ran make for EPD" << endl;
      }
    }

    const Int_t nTracks = mPicoDst->numberOfTracks();
    //    if(cent9 < 0) cout << cent9 << endl;
    //const Double_t reweight = mVecMesonCut->getEventWeight(cent9, refMultCorr);
    //const Int_t nToFMatched = mVecMesonCut->getMatchedToF();
    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      //if(mMode == 4)
      //{
      //  Float_t eta = track->pMom().PseudoRapidity();
      //  if(fabs(eta) < vmsa::mEtaMax && track->gDCA(vx,vy,vz) < 3.0)
      //  {
      //    Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track, mPicoDst);
      //    Float_t dEdx = track->dEdx();
      //    Float_t p = track->pMom().Mag();
      //    mVecMesonHistoManger->FillQA_Detector(dEdx,Mass2,p);
      //  }
      //}
      if(mMode == 4 || mMode == 6)  // fill re-center parameter
      {	
        if(mVecMesonCut->passTrackEP(track,mPicoEvent)) // track cut
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

    //if(mMode == 0) // fill raw EP
    //{
    //  mVecMesonHistoManger->FillQA_Event(vz,refMult);
    //  if(mVecMesonCorrection->passTrackEtaNumRawCut())
    //  {
    //    TVector2 Q2East = mVecMesonCorrection->getQVectorRaw(0); // 0 = eta_gap, 1 = east/west
    //    Float_t Psi2_East = 0.5*TMath::ATan2(Q2East.Y(),Q2East.X());
    //    TVector2 Q2West = mVecMesonCorrection->getQVectorRaw(1); // 0 = eta_gap, 1 = east/west
    //    Float_t Psi2_West = 0.5*TMath::ATan2(Q2West.Y(),Q2West.X());
    //    mVecMesonHistoManger->FillEP_Eta(Psi2_East,Psi2_West);
    //  }
    //  if(mVecMesonCorrection->passTrackFullNumRawCut())
    //  {
    //    TVector2 Q2Full = mVecMesonCorrection->getQVectorRaw(2);
    //    Float_t Psi2_Full = 0.5*TMath::ATan2(Q2Full.Y(),Q2Full.X());
    //    mVecMesonHistoManger->FillEP_Full(Psi2_Full);
    //  }
    //}

    //if(mMode == 1)
    //{
    //  // full event shift parameter
    //  if(mVecMesonCorrection->passTrackFullNumCut())
    //  {
    //    for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
    //    {
    //      TVector2 Psi2Vector_Full_EP = mVecMesonCorrection->calPsi2_Full_EP(k);
    //      mVecMesonProManger->FillEventFull_EP(Psi2Vector_Full_EP,cent9,runIndex,vz_sign,k);
    //    }
    //  }

    //  // eta sub event shift parameter
    //  if(mVecMesonCorrection->passTrackEtaNumCut())
    //  {
    //    for(Int_t k = 0; k < 5; k++)
    //    {
    //      TVector2 Psi2Vector_East_EP = mVecMesonCorrection->calPsi2_East_EP(k);
    //      mVecMesonProManger->FillEventEast_EP(Psi2Vector_East_EP,cent9,runIndex,vz_sign,k);

    //      TVector2 Psi2Vector_West_EP = mVecMesonCorrection->calPsi2_West_EP(k);
    //      mVecMesonProManger->FillEventWest_EP(Psi2Vector_West_EP,cent9,runIndex,vz_sign,k);
    //    }
    //  }
    //}

    //if(mMode == 2) // calculate resolution for eta_sub and random sub event plane
    //{
    //  // calculate Q vector after recentering for Random Sub Event
    //  //cout << "Got to the make section" << endl;
    //  //cout << mUsedTrackCounter << endl;
    //  Int_t iTrack[mUsedTrackCounter];
    //  Float_t ranCounter = (Float_t)mUsedTrackCounter/2.0 - 1;
    //  //cout << "Loop over mUsedTrackCounter" << endl;
    //  for(Int_t i = 0; i < mUsedTrackCounter; i++)
    //  {
    //    iTrack[i] = i;
    //  }
    //  //cout << "Randomly shuffle" << endl;
    //  std::srand(time(0));
    //  std::random_shuffle(iTrack,iTrack+mUsedTrackCounter);

    //  mUsedTrackCounter = 0;
    //  //cout << "Loop over tracks" << endl;
    //  for(Int_t i = 0; i < nTracks; i++) // track loop
    //  {
    //    StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
    //    if(mVecMesonCut->passTrackEP(track,mPicoEvent)) // track cut
    //    {
    //      //cout << "passtrackEP" << endl;
    //      if(mVecMesonCorrection->passTrackFull(track))
    //      {
    //        //cout << "pass track full" << endl;
    //        if((Float_t)iTrack[mUsedTrackCounter] > ranCounter) // Sub Event A
    //        {
    //          //cout << "Sub event A" << endl;
    //          mVecMesonCorrection->addTrack_A(track,cent9,runIndex,vz_sign);
    //        }
    //        else // Sub Event B
    //        {
    //          //cout << "Sub event B" << endl;
    //          mVecMesonCorrection->addTrack_B(track,cent9,runIndex,vz_sign);
    //        }
    //        mUsedTrackCounter++;
    //      }
    //    }
    //  }
    //  mUsedTrackCounter = 0;
    //  //cout << "Calculte resolution" << endl;
    //  // calculate resolution
    //  TVector2 QVecEast = mVecMesonCorrection->getQVector(0);
    //  Float_t Psi2East_ReCenter = 0.5*TMath::ATan2(QVecEast.Y(),QVecEast.X());
    //  Float_t Psi2East_Shift = mVecMesonCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign);

    //  TVector2 QVecWest = mVecMesonCorrection->getQVector(1);
    //  Float_t Psi2West_ReCenter = 0.5*TMath::ATan2(QVecWest.Y(),QVecWest.X());
    //  Float_t Psi2West_Shift = mVecMesonCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign);

    //  TVector2 QVecFull = mVecMesonCorrection->getQVector(2);
    //  Float_t Psi2Full_ReCenter = 0.5*TMath::ATan2(QVecFull.Y(),QVecFull.X());
    //  Float_t Psi2Full_Shift = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);

    //  TVector2 QVecRanA = mVecMesonCorrection->getQVector(3);
    //  Float_t Psi2RanA_ReCenter = 0.5*TMath::ATan2(QVecRanA.Y(),QVecRanA.X());
    //  Float_t Psi2RanA_Shift = mVecMesonCorrection->calShiftAngle2A_EP(runIndex,cent9,vz_sign);

    //  TVector2 QVecRanB = mVecMesonCorrection->getQVector(4);
    //  Float_t Psi2RanB_ReCenter = 0.5*TMath::ATan2(QVecRanB.Y(),QVecRanB.X());
    //  Float_t Psi2RanB_Shift = mVecMesonCorrection->calShiftAngle2B_EP(runIndex,cent9,vz_sign);

    //  if(mVecMesonCorrection->passTrackEtaNumCut())
    //  {
    //    //cout << "pass eta num cut" << endl;
    //    mVecMesonHistoManger->FillEP_Sub(Psi2East_ReCenter,Psi2East_Shift,Psi2West_ReCenter,Psi2West_Shift);
    //    mVecMesonProManger->FillRes_Sub(cent9,Psi2East_Shift,Psi2West_Shift);
    //  }

    //  if(mVecMesonCorrection->passTrackFullNumCut())
    //  {
    //    //cout << "pass track full num cut" << endl;
    //    mVecMesonHistoManger->FillEP_Ran(Psi2RanA_ReCenter,Psi2RanA_Shift,Psi2RanB_ReCenter,Psi2RanB_Shift,Psi2Full_ReCenter,Psi2Full_Shift);
    //    mVecMesonProManger->FillRes_Ran(cent9,Psi2RanA_Shift,Psi2RanB_Shift);
    //  }
    //}

    if(mMode == 4 || mMode == 6)
    { // phi meson
      //if(mVecMesonCorrection->passTrackFullNumCut())
      //{
	// get QVector of sub event
	////TVector2 Q1East = result.EastPhiWeightedQ(1);
        ////TVector2 Q1West = result.WestPhiWeightedQ(1);
	//TVector2 Q2Full = mVecMesonCorrection->getQVector(2); // full 
	//Int_t NumTrackEast = mVecMesonCorrection->getNumTrack(0);
	//Int_t NumTrackWest = mVecMesonCorrection->getNumTrack(1);
	//Int_t NumTrackFull = mVecMesonCorrection->getNumTrack(2);
	//Int_t NumTrackFullEast = mVecMesonCorrection->getNumTrack(3);
	//Int_t NumTrackFullWest = mVecMesonCorrection->getNumTrack(4);

	Float_t Psi1 = result.FullPhiWeightedAndShiftedPsi(1);
        if(mVecMesonCorrection->passTrackFullNumCut())
        {
          Float_t Psi2 = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);
          mVecMesonProManger->FillD12(cent9,Psi1,Psi2); 
        }

	// get N_prim, N_non_prim, N_Tof_match
	////Int_t N_prim = mVecMesonCut->getNpirm();
	////Int_t N_non_prim = mVecMesonCut->getNnonprim();
	////Int_t N_Tof_match = mVecMesonCut->getMatchedToF();

	// pass the event information to StVecMesonTree
	////mVecMesonTree->clearEvent();
	////mVecMesonTree->passEvent(N_prim, N_non_prim, N_Tof_match);
	// pass re-centered event plane to StVecMesonTree
	////mVecMesonTree->passEventPlane(Q1East,Q1West);

	// pass NumOfTrack to StVecMesonTree
	//mVecMesonTree->passNumTrack(NumTrackEast,NumTrackWest,NumTrackFull,NumTrackFullEast,NumTrackFullWest);
        ////if(mFlag_PID == 0)
        ////{
	////  mVecMesonTree->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi1);
        ////}
        ////if(mFlag_PID == 1)
        ////{
	////  mVecMesonTree->MixEvent_Rho(mFlag_ME,mPicoDst,cent9,vz,Psi1);
        ////}
        ////if(mFlag_PID == 2)
        ////{
	////  mVecMesonTree->MixEvent_KStar(mFlag_ME,mPicoDst,cent9,vz,Psi1);
        ////}
      //}
    }
    mVecMesonCorrection->clear();
  }

  return kStOK;
}

