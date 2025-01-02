//Xianglei Zhu 
//Skeleton embedding PicoDst analysis macro with StPicoDstMaker 
//Run it with the wrapper in ACLIC mode, CINT mode for debug ONLY

#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeHelper.h"
#include "TNtuple.h"
#include "TDatime.h"
#include "StarRoot/TUnixTime.h"
#include "StChain.h"
#include "StMessMgr.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StBTofHeader.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/StVecMesonMaker/StUtility.h"
#include "StRoot/StEffKaon/StEffHistManger.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include <map>
//#include <pair>
#include <algorithm>
#endif


void makePicoDstQA_Phi(TString InputFileList, Int_t nEvents = 0, TString OutputFile = "test.histo", TString jobId = "1", Int_t mEnergy = 4, Int_t mGid = 11, Int_t mPid = 0);

void makePicoDstQA_Phi(TString InputFileList, Int_t nEvents, TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid, Int_t mPid) 
{
 
  // Load libraries for CINT mode
#ifdef __CINT__
  gROOT  -> Macro("loadMuDst.C");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StVecMesonMaker");
  gSystem->Load("StEffKaon");
  gSystem->Load("StUtility");
#endif

  // List of member links in the chain
  StChain*                    chain  =  new StChain ;

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead,InputFileList,"picoDst");

  StUtility* mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index 
  //mUtility->initEventPlane(); // initialize std::map for event planes 

  StRefMultCorr *mRefMultCorr = new StRefMultCorr("refmult");

  StVecMesonCut* mVecMesonCut = new StVecMesonCut(mEnergy);

  StEffHistManger* mEffHistManger = new StEffHistManger();

  if ( nEvents == 0 )  nEvents = 1000000000 ;       // Take all events in nFiles if nEvents = 0

  // ---------------- modify here according to your QA purpose --------------------------
  TFile *tags_output = new TFile( OutputFile+"_"+jobId+".root" , "recreate" ) ;
  tags_output->cd();

  // ---------------- end of histogram and tree booking --------------------------------

  // Loop over the links in the chain
  Int_t iInit = chain -> Init() ;
  if (iInit) chain->FatalErr(iInit,"on init");

  mEffHistManger->InitHist();
  
  // chain -> EventLoop(1,nEvents) ;  //will output lots of useless debugging info.
  Int_t istat = 0, i = 1;
  while (i <= nEvents && istat != 2) {
     if(i%1000==0)cout << endl << "== Event " << i << " start ==" << endl;
     chain->Clear();
     //cout << "Event = " << i << endl;
     istat = chain->Make(i);

     if (istat == 2)
	  cout << "Last  event processed. Status = " << istat << endl;
     if (istat == 3)
	  cout << "Error event processed. Status = " << istat << endl;
     i++;

     if(istat != kStOK)continue; //skip those suspectible events
     
  // ---------------- modify here according to your QA purpose --------------------------
     //cout<<"In event #. "<<i-1<<" Maker status "<<istat<<endl;

     StPicoDst* mPicoDst = picoDstMaker->picoDst();
     if(!mPicoDst) {
	  LOG_WARN << " No PicoDst " << endm; continue;
     }

     StPicoEvent* mPicoEvent = mPicoDst->event();
     if(!mPicoEvent) {
	  LOG_WARN << " No PicoEvent " << endm; continue;
     }

     //triggerdd
     //if ( ! mPicoEvent->isTrigger(640001) && ! mPicoEvent->isTrigger(640011) && ! mPicoEvent->isTrigger(640021) && ! mPicoEvent->isTrigger(640031) && ! mPicoEvent->isTrigger(640041) && ! mPicoEvent->isTrigger(640051) ) continue ; // BES-II 19.6 GeV
    
     Int_t runId = mPicoEvent->runId();
     Int_t eventId = mPicoEvent->eventId();
     Int_t refMult = mPicoEvent->refMult();
     Float_t vx = mPicoEvent->primaryVertex().x();
     Float_t vy = mPicoEvent->primaryVertex().y();
     Float_t vz = mPicoEvent->primaryVertex().z();
     Float_t zdcX = mPicoEvent->ZDCx();
     const unsigned short nBTofMatched = mPicoEvent->nBTOFMatch();

     mRefMultCorr->init(runId);

     if(mRefMultCorr->isBadRun( runId ))
     {
       LOG_ERROR << "Bad Run! Skip!" << endm; continue;
     }
   
     bool isPileUpEvent = false;
     // IMPORTANT: vertex position is needed for Au+Au 19.6 GeV 2019
     if (mRefMultCorr->isPileUpEvent( refMult, nBTofMatched, vz ) ) isPileUpEvent = true;
     mRefMultCorr->initEvent(refMult,vz,zdcX);

     const Int_t cent9 = mRefMultCorr->getCentralityBin9();       // 0: 70-80%, 1: 60-70%, ..., 6: 10-20%, 7: 5-10%, 8: 0-5%
 
     if(cent9 < -0.5) continue;
     if(isPileUpEvent) continue;
     if(!mVecMesonCut->passEventCut(mPicoDst)) continue;

     const int runIndex = mUtility->findRunIndex(runId);
      
     if(runIndex < 0)
     {
       LOG_ERROR << "Could not find this run Index from StUtility! Skip!" << endm; continue;
     }

     const float ep_west = 0.0;
     const float ep_east = 0.0;
     const float ep_full = 0.0;
     //const float ep_west = mUtility->findEPwest(eventId);
     //const float ep_east = mUtility->findEPeast(eventId);
     //const float ep_full = mUtility->findEPfull(eventId);
     //if( ep_west < -900.0 || ep_east < -900.0 || ep_full < -900.0 )
     //{
     //  LOG_ERROR << "Could not find this event plane from StUtility! Skip!" << endm; continue;
     //}

     //Vz
     //if ( fabs(mPicoEvent->primaryVertex().Z()) > 70.0 ) continue ;
     //Vr
     //if ( mPicoEvent->primaryVertex().Perp() > 2.0 ) continue ;
     
     //fill MC histograms
     //The MC arrays in PicoDst
     Int_t NoMuMcVertices = mPicoDst->numberOfMcVertices();
     Int_t NoMuMcTracks = mPicoDst->numberOfMcTracks();
     //LOG_INFO <<"# of MC tracks = "<< NoMuMcTracks <<" # of MC vertices = "<< NoMuMcVertices << endm;
     if (! NoMuMcVertices || ! NoMuMcTracks) {
	  //LOG_WARN << "Ev. " << i  << " has no MC information ==> skip it" << endm;
	  continue;
     }
     Int_t nMc = 0;

     // Loop for MC tracks
     /*for(Int_t itrk=0; itrk<NoMuMcTracks; itrk++){
	  StPicoMcTrack *mcTrack = (StPicoMcTrack *) mPicoDst->mcTrack(itrk);
	  if (! mcTrack) continue;

	  // Select only Triggered Mc Vertex, i.e. the MC track should originate from PV (IdVx=1)
	  Int_t IdVx = mcTrack->idVtxStart();
	  if (IdVx != 1) continue;

	  const int Gid = mcTrack->geantId();
          //nMc++;  // # of MC tracks
	  //if(Gid==11){//k+
	  if(Gid==mGid){//k-
		//hPtMc->Fill(mcTrack->p().Perp());
		//hPhiMc->Fill(mcTrack->p().Phi());
		//hEtaMc->Fill(mcTrack->p().PseudoRapidity());
		//if(fabs(mcTrack->p().PseudoRapidity())<1.5)hSelPtMc->Fill(mcTrack->p().Perp()); //This simply limits the eta range for selected tracks
	  }
	  else {
	     LOG_WARN << "Gid: "<<Gid<<" in Ev. "<<i<<endm;
             continue;
	  }
         
          double P   = mcTrack->p().Mag();
          double Pt  = mcTrack->p().Perp();
          double Phi = mcTrack->p().Phi();
          double Eta = mcTrack->p().PseudoRapidity();

          //-------------------------McKaon-----------------------------------------------------
          if( Eta > vmsa::mEtaMax ) continue; // eta cut
          if(!(Pt > vmsa::mGlobPtMin && P < vmsa::mPrimMomMax) ) continue; // eta cut

          mEffHistManger->FillHistMc(cent9,Pt,Eta,Phi,0.0,0.0);
          //-------------------------McKaon----------------------------------------------------- 
     }*/

     //fill Event QA histograms
     Int_t nTracks = mPicoDst->numberOfTracks();
     StPicoTrack* ptrack ; 

     std::map< int, std::pair<StPicoMcTrack*, StPicoMcTrack*> > phiMesonMcDaughters;
     std::map< int, std::pair<StPicoTrack*, StPicoTrack*> > phiMesonRcDaughters;
 
     for(Int_t i=0; i<nTracks; i++)                // Main loop for Iterating over tracks
     {
	  ptrack = mPicoDst->track(i);  // Pointer to a track
	  if(!ptrack->isPrimary()) continue;

	  if (ptrack->idTruth() <= 0 || ptrack->idTruth() > NoMuMcTracks) {
	     //cout << "Illegal idTruth " << ptrack->idTruth() << " The track is ignored" << endl;
	     continue;
	  }
	  StPicoMcTrack *mcTrack = (StPicoMcTrack *) mPicoDst->mcTrack(ptrack->idTruth()-1);
	  if (!mcTrack) {
	     LOG_WARN << "Inconsistency in mcArray(1), ignored" << endm;
	     continue;
	  }
	  if (mcTrack->id() != ptrack->idTruth()) {
	     LOG_WARN << "Mismatched idTruth " << ptrack->idTruth() << " and mcTrack Id " <<  mcTrack->id() 
		  << " this track is ignored" <<  endm;
	  }
	  Int_t idMcVx = mcTrack->idVtxStart();
          Int_t parentGeantId = -1;
          Int_t idMcTrack = -1;
	  while (idMcVx != 1) {
	     StPicoMcVertex *mcVertex = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVx-1);
	     idMcTrack = mcVertex->idOfParentTrack();
             //cout << "idMcTrack = " << idMcTrack << endl;
	     if (! idMcTrack) break;
	     StPicoMcTrack *mcTrackP = (StPicoMcTrack *) mPicoDst->mcTrack(idMcTrack-1);
             parentGeantId = mcTrackP->geantId();
	     idMcVx = mcTrackP->idVtxStart();
	     if (! idMcVx) break;
	  }
	  if (idMcVx != 1) continue; //this MC track is not eventually originated from PV

          if (parentGeantId != 10151) continue;
          if (idMcTrack == -1) continue;

	  //if(mcTrack->geantId() != 11) continue; //geantId cut for Kplus = 11 and Kminus = 12
	  if(mcTrack->geantId() != 11 && mcTrack->geantId() != 12) continue; //geantId cut for Kplus = 11 and Kminus = 12
	  //if(mcTrack->idVtxStart() != 1) {
	  //   LOG_WARN<<"mc track may not directly originate from PV!"<<endm;
	  //}

          double P   = mcTrack->p().Mag();
          double Pt  = mcTrack->p().Perp();
          //double Phi = mcTrack->p().Phi();
          double Eta = mcTrack->p().PseudoRapidity();
          //double Y = mcTrack->p().Rapidity();
   
          //double gRcPt = ptrack->gMom().Perp();
          //double pRcPt = ptrack->pMom().Perp();
 
          if(ptrack->qaTruth()<50.) continue;

          //cout << "nsigma_kaon = " << ptrack->nSigmaKaon() << endl;    

          //-------------------------McKaon-----------------------------------------------------
          if( fabs(Eta) > vmsa::mEtaMax ) continue; // eta cut
          //if(!(Pt > vmsa::mGlobPtMin && P < vmsa::mPrimMomMax) ) continue; // pt cut

          if(mcTrack->geantId() == 11) 
          {
            phiMesonMcDaughters[idMcTrack].first = mcTrack;
            phiMesonRcDaughters[idMcTrack].first = ptrack;
            //cout << "index = " << idMcTrack << " kplus mc" << endl;
          }
          if(mcTrack->geantId() == 12) 
          {
            phiMesonMcDaughters[idMcTrack].second = mcTrack;
            phiMesonRcDaughters[idMcTrack].second = ptrack;
            //cout << "index = " << idMcTrack << " kminus mc" << endl;
          }

          //mEffHistManger->FillHistMc(cent9,Pt,Eta,Phi,ep_west,ep_east);
          //-------------------------McKaon----------------------------------------------------- 

          //if( !mVecMesonCut->passTrackMeson(ptrack, mPicoEvent, mPid) ) continue;
          //mEffHistManger->FillHistRc(cent9,Pt,Eta,Phi,ep_west,ep_east);
          //mEffHistManger->FillHistPt(cent9,Pt,gRcPt,pRcPt);   
          //cout << "kaon m^2 = " << mVecMesonCut->getPrimaryMass2(ptrack, mPicoDst) << endl;    

	  //if(ptrack->qaTruth()<50.) continue; // Quality of MC track (probably out of 100)

	  //if(ptrack->nHits()<=15)continue; // nHits cut
	  //if(ptrack->flag()<=0)continue;
	  //if(abs(ptrack->charge())!=1) continue; //if the charge of the track is not equal to 1

	  //TVector3 p = ptrack->pMom();
	  //hPhi->Fill(p.Phi());
	  //hPt->Fill(p.Perp());
	  //hEta->Fill(p.PseudoRapidity());
	  //if(fabs(p.PseudoRapidity())<1.5)hSelPt->Fill(p.Perp());
	  //end of the filling  
     }
   
     //std::sort( phiIndex.begin(), phiIndex.end() );
     //phiIndex.erase( std::unique( phiIndex.begin(), phiIndex.end() ), phiIndex.end() );
 
     for(std::map< int, std::pair<StPicoMcTrack*,StPicoMcTrack*> >::iterator phi = phiMesonMcDaughters.begin(); phi != phiMesonMcDaughters.end(); phi++) 
     {
       int index = (phi->first);
       StPicoMcTrack *mcTrackKP = phiMesonMcDaughters[index].first;
       StPicoMcTrack *mcTrackKM = phiMesonMcDaughters[index].second;
       if(mcTrackKP == nullptr || mcTrackKM == nullptr) continue;
       double KP_P   = mcTrackKP->p().Mag();
       double KP_Pt  = mcTrackKP->p().Perp();
       double KP_Phi = mcTrackKP->p().Phi();
       double KP_Eta = mcTrackKP->p().PseudoRapidity();
       double KM_P   = mcTrackKM->p().Mag();
       double KM_Pt  = mcTrackKM->p().Perp();
       double KM_Phi = mcTrackKM->p().Phi();
       double KM_Eta = mcTrackKM->p().PseudoRapidity();


       TLorentzVector l_kp, l_km;
       l_kp.SetPtEtaPhiM(KP_Pt,KP_Eta,KP_Phi,vmsa::mMassKaon);
       l_km.SetPtEtaPhiM(KM_Pt,KM_Eta,KM_Phi,vmsa::mMassKaon);
       TLorentzVector lphi = l_kp + l_km;        

       //cout << "index = " << index << " mc phi candidate" << endl;
    
       if(TMath::Abs(lphi.Rapidity()) > 1.0) continue;

       //cout << "index = " << index << " mc phi candidate passed rapidity cut" << endl;

       mEffHistManger->FillHistMc(cent9,lphi.Pt(),lphi.Eta(),lphi.Phi(),ep_west,ep_east);

       StPicoTrack *ptrackKP = phiMesonRcDaughters[index].first;
       StPicoTrack *ptrackKM = phiMesonRcDaughters[index].second;
       if(ptrackKP == nullptr || ptrackKM == nullptr) continue;
       double rKP_P   = ptrackKP->pMom().Mag();
       double rKP_Pt  = ptrackKP->pMom().Perp();
       double rKP_Phi = ptrackKP->pMom().Phi();
       double rKP_Eta = ptrackKP->pMom().PseudoRapidity();
       double rKM_P   = ptrackKM->pMom().Mag();
       double rKM_Pt  = ptrackKM->pMom().Perp();
       double rKM_Phi = ptrackKM->pMom().Phi();
       double rKM_Eta = ptrackKM->pMom().PseudoRapidity();

       //cout << "index = " << index << " rc phi candidate" << endl;

       if( !mVecMesonCut->passTrackMeson(ptrackKP, mPicoEvent, mPid) ) continue;
       if( !mVecMesonCut->passTrackMeson(ptrackKM, mPicoEvent, mPid) ) continue;

       //cout << "index = " << index << " rc phi pass" << endl;

       TLorentzVector rl_kp, rl_km;
       rl_kp.SetPtEtaPhiM(rKP_Pt,rKP_Eta,rKP_Phi,vmsa::mMassKaon);
       rl_km.SetPtEtaPhiM(rKM_Pt,rKM_Eta,rKM_Phi,vmsa::mMassKaon);
       TLorentzVector rlphi = rl_kp + rl_km;

       if(TMath::Abs(rlphi.Rapidity()) > 1.0) continue;
       //cout << "index = " << index << " mc phi candidate passed rapidity cut" << endl;

       mEffHistManger->FillHistRc(cent9,lphi.Pt(),lphi.Eta(),lphi.Phi(),ep_west,ep_east);

       //cout << "azimuthal angle = " << lphi.Phi() << endl;
     }

     phiMesonMcDaughters.clear();
     phiMesonRcDaughters.clear();
  }
  //mEffHistManger->CalEfficiency();
  //mEffHistManger->CalEffPtEtaPhi();

  tags_output->cd();
  mEffHistManger->WriteHist();

  //hEffPt->Divide(hSelPt,hSelPtMc,1,1,"B");

  if (nEvents > 1) chain -> Finish() ;

  if(tags_output!=NULL) tags_output -> Write() ;
  if(tags_output!=NULL) tags_output -> Close() ;
  //flush(tags_output);
  delete tags_output;

  // Cleanup
  delete chain ;
}
