#include "StRoot/StVecMesonMaker/StVecMesonTree.h"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StMesonEvent/StMesonEvent.h"
#include <vector>
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "StarClassLibrary/StThreeVector.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "TMath.h"
#include "TObject.h"
#include "TVector3.h"
#include "TRandom3.h"

ClassImp(StVecMesonTree)

//------------------------------------------------------------------------------------------------------------------
StVecMesonTree::StVecMesonTree(Int_t energy)
{
  mEnergy = energy;
}

StVecMesonTree::~StVecMesonTree()
{
  /* */
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::InitPhi()
{
  
  gRandom->SetSeed();
 
  cout << "StVecMesonTree::InitPhi()" << endl;
  cout << "StVecMesonCut" << endl;
  mVecMesonCut = new StVecMesonCut(mEnergy);
  cout << "Creating Histograms" << endl;
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.98,1.08);
  HistName = "K_dEdx_Rig";
  h_KdEdxRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,400,0,40);
  HistName = "K_Mass2_Rig";
  h_KM2Rig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,100,-0.25,0.5);
  HistName = "K_InvBeta_Rig";
  h_KInvBetaRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,300,-0.0,3.0);

  cout << "clearing everything" << endl;
  for(Int_t cent = 0; cent < vmsa::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < vmsa::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < vmsa::Bin_Phi_Psi; phi_psi++)
      {
        cout << "cent = " << cent << ",  vz = " << vz << ",  phi_psi = " << phi_psi << endl;
        cout << "Set event counter to 0 " << endl;
        mEventCounter2[cent][vz][phi_psi] = 0;
        cout << "clear_phi" << endl;
	clear_phi(cent,vz,phi_psi);

        cout << "Creating 2 StMesonEvents" << endl;
        mMesonEvents[cent][vz][phi_psi][0] = new StMesonEvent();
        mMesonEvents[cent][vz][phi_psi][1] = new StMesonEvent();

        mMixed[cent][vz][phi_psi][0] = 0;
        mMixed[cent][vz][phi_psi][1] = 0;
      }
    }
  }
  cout << "Creating a TTree" << endl;
  mTree = new TTree(vmsa::vm_tree[0].Data(),vmsa::vm_tree[0].Data());
  cout << "Creating a branch" << endl;
  mMesonEvent = new StMesonEvent();
  mTree->Branch(vmsa::vm_branch[0].Data(),"StMesonEvent",&mMesonEvent);
  cout << "SetAutoSave" << endl;
  //mTree->SetAutoSave(5000000);
}


void StVecMesonTree::WriteMass2Phi()
{
  h_Mass2->Write();
  h_KdEdxRig->Write();
  h_KM2Rig->Write();
  h_KInvBetaRig->Write(); 
  mTree->Write("",TObject::kOverwrite);
}

void StVecMesonTree::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  mMixed[cent9][Bin_vz][Bin_Psi2][0] = 0;
  mMixed[cent9][Bin_vz][Bin_Psi2][1] = 0;
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunId[cent9][Bin_vz][Bin_Psi2].clear();
  mEventId[cent9][Bin_vz][Bin_Psi2].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2East[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2West[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2Full[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackFull[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t Bin_Event = 0; Bin_Event <= vmsa::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
      mHelix[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigma[key].clear();
      mEventBin[key].clear();
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StVecMesonTree::size_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  LOG_INFO << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << Bin_vz << ", Psi2 = " << Bin_Psi2 << endm;
  LOG_INFO << "Buffer_depth = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endm;

  LOG_INFO << "Size of primaryVertex = " << mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].size() << endm;;
  LOG_INFO << "Size of refMult       = " << mRefMult[cent9][Bin_vz][Bin_Psi2].size() << endm;;
  LOG_INFO << "---------------------------------------------------------------------------" << endm;

  for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
  {
    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
    LOG_INFO << "Event Number " << Bin_Event << ":" << endm; 
    LOG_INFO << "Positive Kaon:" << endm;
    LOG_INFO << "  Size of Helix        = " << mHelix[key].size() << endm;;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigma       = " << mNSigma[key].size() << endm;
    LOG_INFO << "  Size of Charge       = " << mCharge[key].size() << endm;

    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
    LOG_INFO << "Negative Kaon:" << endm;
    LOG_INFO << "  Size of Helix        = " << mHelix[key].size() << endm;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigma       = " << mNSigma[key].size() << endm;
    LOG_INFO << "  Size of Charge       = " << mCharge[key].size() << endm;
    LOG_INFO << "---------------------------------------------------------------------------" << endm;
  }
}


void StVecMesonTree::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2, Int_t &EventToMix) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 1) // mixed event
  {
    // Flip the coin for which event to select
    double randomValue = gRandom->Uniform();
    EventToMix = (randomValue < 0.5) ? 0 : 1; // Event 0 or 1 
    //cout << "EventToMix = " << EventToMix << endl;
    
    MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,EventToMix,0);
    MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,EventToMix,1);
    
    MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,2,0); // 3rd event, repeatedly replaced by new events  
    MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,2,1); // 3rd event, repeatedly replaced by new events
    
    Int_t Bin_Event = EventToMix;
    
    //Don't need to clear this anymore because we are going to save all of these guys
    //mMesonEvent->clearTrackList();

    //cout << "What event are we mixing? " << Bin_Event << endl;
   
    if(mMixed[cent9][Bin_vz][Bin_Psi2][Bin_Event] == 0)
    {
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      
      // QVector
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      
      // Number of Tracks
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      
      mMixed[cent9][Bin_vz][Bin_Psi2][Bin_Event] = 1;
      //cout << "Set mMixed to 1" << endl;
    }

    TLorentzVector ltrackA, ltrackB;
    // start to mix events
    // mix K+ candidates from A event with K- candidates from B event
    for(Int_t n_kplus = 0; n_kplus < mHelix[key_A_plus].size(); n_kplus++) // first track loop over K+ candidates from event selected to mix
    {
      StThreeVectorF p_vecA(mHelix[key_A_plus][n_kplus].cat(mHelix[key_A_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])));  // primary momentum
      p_vecA *= mMomentum[key_A_plus][n_kplus];
      ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+
    
      for(Int_t n_kminus = 0; n_kminus < mHelix[key_B_minus].size(); n_kminus++) // second track loop over K- candidates from the 3rd event 
      {
        StThreeVectorF p_vecB(mHelix[key_B_minus][n_kminus].cat(mHelix[key_B_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][2])));  // primary momentum
        p_vecB *= mMomentum[key_B_minus][n_kminus];
        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-
    
        TLorentzVector trackAB      = ltrackA+ltrackB;
        Double_t InvMassAB          = trackAB.M();
        Double_t pt = trackAB.Perp();
        // fill phi candidate background into mTree
        if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
        {
          mMesonTrack = mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->createTrack();
          mMesonTrack->setMass2A(mMass2[key_A_plus][n_kplus]); // K+
          mMesonTrack->setMass2B(mMass2[key_B_minus][n_kminus]); // K-
          mMesonTrack->setNSigA(mNSigma[key_A_plus][n_kplus]); // K+
          mMesonTrack->setNSigB(mNSigma[key_B_minus][n_kminus]); // K-
          mMesonTrack->setDcaA(mDca[key_A_plus][n_kplus]); // K+
          mMesonTrack->setDcaB(mDca[key_B_minus][n_kminus]); // K-
          mMesonTrack->setTrackA(ltrackA); // K+
          mMesonTrack->setTrackB(ltrackB); // K-
          mMesonTrack->setFlagA(mEventBin[key_A_plus][n_kplus]); // K+
          mMesonTrack->setFlagB(2); // K-
          h_Mass2->Fill(pt,InvMassAB);
        }
      }
    }
    
    // mix K- candidates from A event with K+ candidates from B event
    for(Int_t n_kminus = 0; n_kminus < mHelix[key_A_minus].size(); n_kminus++) // first track loop over K- candidates from event selected to mix
    {
      StThreeVectorF p_vecA(mHelix[key_A_minus][n_kminus].cat(mHelix[key_A_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])));  // primary momentum
      p_vecA *= mMomentum[key_A_minus][n_kminus];
      ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-
    
      for(Int_t n_kplus = 0; n_kplus < mHelix[key_B_plus].size(); n_kplus++) // second track loop over K+ candidates from the 3rd event
      {
        StThreeVectorF p_vecB(mHelix[key_B_plus][n_kplus].cat(mHelix[key_B_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][2])));  // primary momentum
        p_vecB *= mMomentum[key_B_plus][n_kplus];
        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+
    
        TLorentzVector trackAB      = ltrackA+ltrackB;
        Double_t InvMassAB          = trackAB.M();
        Double_t pt = trackAB.Perp();
    
        // fill phi candidate background into mTree
        if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
        {
          mMesonTrack = mMesonEvents[cent9][Bin_vz][Bin_Psi2][Bin_Event]->createTrack();
          mMesonTrack->setMass2A(mMass2[key_B_plus][n_kplus]); // K+
          mMesonTrack->setMass2B(mMass2[key_A_minus][n_kminus]); // K-
          mMesonTrack->setNSigA(mNSigma[key_B_plus][n_kplus]); // K+
          mMesonTrack->setNSigB(mNSigma[key_A_minus][n_kminus]); // K-
          mMesonTrack->setDcaA(mDca[key_B_plus][n_kplus]); // K+
          mMesonTrack->setDcaB(mDca[key_A_minus][n_kminus]); // K-
          mMesonTrack->setTrackA(ltrackB); // K+
          mMesonTrack->setTrackB(ltrackA); // K-
          mMesonTrack->setFlagA(2); // K+
          mMesonTrack->setFlagB(mEventBin[key_A_minus][n_kminus]); // K-
          h_Mass2->Fill(pt,InvMassAB);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::MixEvent_Phi(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2)
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz, Bin_Psi2;

  Float_t vz_start = vmsa::mVzMaxMap[mEnergy];
  Float_t vz_bin = 2*vz_start/vmsa::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/vmsa::Bin_Phi_Psi;

  for(Int_t i = 0; i < vmsa::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < vmsa::Bin_Phi_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  //cout << "Bin Event = " << Bin_Event << endl;

  const Double_t MAGFIELDFACTOR = kilogauss;
  //cout << "Tesla     = " << tesla << endl;
  //cout << "KiloGauss = " << kilogauss << endl;
  //cout << "Gauss     = " << gauss << endl;
  const Int_t nTracks = pico->numberOfTracks();

  // store Event Information
  StThreeVectorF primVer; 
  primVer.setX(event->primaryVertex().x());
  primVer.setY(event->primaryVertex().y());
  primVer.setZ(event->primaryVertex().z());

  //if(Bin_Event < 2) mMixed[cent9][Bin_vz][Bin_Psi2][0] = 0;
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(primVer));
  mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  mEventId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->eventId()));
  mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  mQ2East[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2East));
  mQ2West[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2West));
  mQ2Full[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2Full));
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaEast));
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaWest));
  mNumTrackFull[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFull));
  mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullEast));
  mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullWest));

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);
    if(mVecMesonCut->passTrackMeson(track,event,0))
    {
      Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track,pico);
      //Float_t scale_nSigma_factor = vmsa::mSigScaleMap[mEnergy];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().Mag();
      Float_t eta = track->pMom().PseudoRapidity();
      Float_t Mass2_low = 0.16;
      Float_t Mass2_up = 0.36;

      Int_t charge = 0; // k+
      if(Polarity < 0) charge = 1; // k-

      if(mVecMesonCut->passSigKaonCut(track,0)) 
      {
	//if(
	//    (momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
	//    || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)
	//  )
        //cout << "fabs(eta) = " << fabs(eta) << endl;
    	if( fabs(eta) <= 1.0 && (Mass2 < Mass2_low || Mass2 > Mass2_up) ) continue;
        if( fabs(eta) > 1.0) continue;
	
        StThreeVectorD primMom; 
        primMom.setX(track->pMom().x());
        primMom.setY(track->pMom().y());
        primMom.setZ(track->pMom().z());
        StThreeVectorD primVer; 
        primVer.setX(event->primaryVertex().x());
        primVer.setY(event->primaryVertex().y());
        primVer.setZ(event->primaryVertex().z());
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	mEventBin[key].push_back(static_cast<Int_t>(Bin_Event)); // mass2
	mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
	mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
	mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaKaon()))); // nSigma
	mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
        //cout << "MagneticField = " << event->bField() << "    MAGFIELDFACTOR = " << MAGFIELDFACTOR << "    product = " << event->bField()*MAGFIELDFACTOR << endl; 
	mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
        h_KdEdxRig->Fill(momentum*Polarity,track->dEdx());
        h_KM2Rig->Fill(momentum*Polarity,Mass2);
        h_KInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
      }
    }
  }
 
  //cout << "Added all of the trakcks" << endl;

  if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] <= 1) mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  //cout << "EventCounter: " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endl;

  //cout << "Flag_ME = " << Flag_ME << endl; 

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] >= vmsa::Buffer_depth)
    {
      
      //cout << "About to doPhi()" << endl;
      Int_t EventToMix = -1;
      
      doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2,EventToMix);
      //cout << "Referenced EventToMix = " << EventToMix << endl << endl;
  
      //mMixed[cent9][Bin_vz][Bin_Psi2].pop_back();
      mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].pop_back();
      mRefMult[cent9][Bin_vz][Bin_Psi2].pop_back();
      mCentrality[cent9][Bin_vz][Bin_Psi2].pop_back();
      mRunId[cent9][Bin_vz][Bin_Psi2].pop_back();
      mEventId[cent9][Bin_vz][Bin_Psi2].pop_back();
      mN_prim[cent9][Bin_vz][Bin_Psi2].pop_back();
      mN_non_prim[cent9][Bin_vz][Bin_Psi2].pop_back();
      mN_Tof_match[cent9][Bin_vz][Bin_Psi2].pop_back();
      mZDCx[cent9][Bin_vz][Bin_Psi2].pop_back();
      mBBCx[cent9][Bin_vz][Bin_Psi2].pop_back();
      mVzVpd[cent9][Bin_vz][Bin_Psi2].pop_back();
      mNumTracks[cent9][Bin_vz][Bin_Psi2].pop_back();
      mQ2East[cent9][Bin_vz][Bin_Psi2].pop_back();
      mQ2West[cent9][Bin_vz][Bin_Psi2].pop_back();
      mQ2Full[cent9][Bin_vz][Bin_Psi2].pop_back();
      mNumTrackEast[cent9][Bin_vz][Bin_Psi2].pop_back();
      mNumTrackWest[cent9][Bin_vz][Bin_Psi2].pop_back();
      mNumTrackFull[cent9][Bin_vz][Bin_Psi2].pop_back();
      mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].pop_back();
      mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].pop_back();

      //mEventCounter2[cent9][Bin_vz][Bin_Psi2]--;
    
      MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,EventToMix,0); // event selected to mix, repeatedly replaced by new events  
      MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,EventToMix,1); // event selected to mix, repeatedly replaced by new events
      MEKey key_plus_3rd  = MEKey(cent9,Bin_vz,Bin_Psi2,2,0); // 3rd event, repeatedly replaced by new events  
      MEKey key_minus_3rd = MEKey(cent9,Bin_vz,Bin_Psi2,2,1); // 3rd event, repeatedly replaced by new events

      for(Int_t n_kminus = 0; n_kminus < mHelix[key_minus_3rd].size(); n_kminus++) // Clear container of K- candidates from the 3rd event 
      {
	mEventBin[key_minus].push_back(mEventBin[key_minus_3rd][n_kminus]);
	mMass2[key_minus].push_back(mMass2[key_minus_3rd][n_kminus]);
	mDca[key_minus].push_back(mDca[key_minus_3rd][n_kminus]);
	mNHitsFit[key_minus].push_back(mNHitsFit[key_minus_3rd][n_kminus]);
	mNSigma[key_minus].push_back(mNSigma[key_minus_3rd][n_kminus]);
	mHelix[key_minus].push_back(mHelix[key_minus_3rd][n_kminus]);
	mMomentum[key_minus].push_back(mMomentum[key_minus_3rd][n_kminus]);

	mEventBin[key_minus_3rd].pop_back();
	mMass2[key_minus_3rd].pop_back();
	mDca[key_minus_3rd].pop_back();
	mNHitsFit[key_minus_3rd].pop_back();
	mNSigma[key_minus_3rd].pop_back();
	mHelix[key_minus_3rd].pop_back();
	mMomentum[key_minus_3rd].pop_back();
      }
      for(Int_t n_kplus = 0; n_kplus < mHelix[key_plus_3rd].size(); n_kplus++) // Clear container of K+ candidates from the 3rd event 
      {
	mEventBin[key_plus].push_back(mEventBin[key_plus_3rd][n_kplus]);
	mMass2[key_plus].push_back(mMass2[key_plus_3rd][n_kplus]);
	mDca[key_plus].push_back(mDca[key_plus_3rd][n_kplus]);
	mNHitsFit[key_plus].push_back(mNHitsFit[key_plus_3rd][n_kplus]);
	mNSigma[key_plus].push_back(mNSigma[key_plus_3rd][n_kplus]);
	mHelix[key_plus].push_back(mHelix[key_plus_3rd][n_kplus]);
	mMomentum[key_plus].push_back(mMomentum[key_plus_3rd][n_kplus]);

	mEventBin[key_plus_3rd].pop_back();
	mMass2[key_plus_3rd].pop_back();
	mDca[key_plus_3rd].pop_back();
	mNHitsFit[key_plus_3rd].pop_back();
	mNSigma[key_plus_3rd].pop_back();
	mHelix[key_plus_3rd].pop_back();
	mMomentum[key_plus_3rd].pop_back();
      }
    }
  }
}

void StVecMesonTree::FillTree_Phi()
{
  cout << "About to fill the TTree" << endl;
  for(Int_t cent = 0; cent < vmsa::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < vmsa::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < vmsa::Bin_Phi_Psi; phi_psi++)
      {
        cout << "cent: " << cent << ", vz: " << vz << ", phi-Psi: " << phi_psi << endl;
        mMesonEvent = mMesonEvents[cent][vz][phi_psi][0];
        //cout << "Assigned mMesonEvent for event 0" << endl;
        cout << "Value of mMixed = " << mMixed[cent][vz][phi_psi][0] << endl;
        if(mMixed[cent][vz][phi_psi][0] != 0) mTree->Fill(); // Fill only if this event has been mixed with other events
        //cout << "Filled TTree with event 0" << endl;
        mMesonEvent = mMesonEvents[cent][vz][phi_psi][1];
        //cout << "Assigned mMesonEvent for event 1" << endl;
        cout << "Value of mMixed = " << mMixed[cent][vz][phi_psi][1] << endl;
        if(mMixed[cent][vz][phi_psi][1] != 0) mTree->Fill(); // Fill only if this event has been mixed with other events
        //cout << "Filled TTree with event 1" << endl;
      }
    }
  }
  cout << "Filled the TTree" << endl;
}


//------------------------------------------------------------------------------------------------------------------

// pass event information from Maker
void StVecMesonTree::clearEvent()
{
  mNumber_prim = 0;
  mNumber_non_prim = 0;
  mNumber_Tof_match = 0;

  mQVector2East.Set(-999.9,-999.9);
  mQVector2West.Set(-999.9,-999.9);
  mQVector2Full.Set(-999.9,-999.9);

  mTrackEtaEast  = 0;
  mTrackEtaWest  = 0;
  mTrackFull     = 0;
  mTrackFullEast = 0;
  mTrackFullWest = 0;
}

void StVecMesonTree::passEvent(Int_t N_prim, Int_t N_non_prim, Int_t N_Tof_match)
{
  mNumber_prim = N_prim;
  mNumber_non_prim = N_non_prim;
  mNumber_Tof_match = N_Tof_match;
}

void StVecMesonTree::passEventPlane(TVector2 Q2East, TVector2 Q2West, TVector2 Q2Full)
{
  mQVector2East = Q2East;
  mQVector2West = Q2West;
  mQVector2Full = Q2Full;
}

void StVecMesonTree::passNumTrack(Int_t NumTrackEast, Int_t NumTrackWest, Int_t NumTrackFull, Int_t NumTrackFullEast, Int_t NumTrackFullWest)
{
  mTrackEtaEast  = NumTrackEast;
  mTrackEtaWest  = NumTrackWest;
  mTrackFull     = NumTrackFull;
  mTrackFullEast = NumTrackFullEast;
  mTrackFullWest = NumTrackFullWest;
}
//------------------------------------------------------------------------------------------------------------------
