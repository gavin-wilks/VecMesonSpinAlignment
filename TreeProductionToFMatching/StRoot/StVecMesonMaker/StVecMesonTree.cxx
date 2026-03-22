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
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
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
        //cout << "cent = " << cent << ",  vz = " << vz << ",  phi_psi = " << phi_psi << endl;
        mEventCounter2[cent][vz][phi_psi] = 0;
	clear_phi(cent,vz,phi_psi);
      }
    }
  }

  mMesonEventSE = new StMesonEvent();
  std::string SEname = Form("%s_SE",vmsa::vm_tree[0].Data());
  //mTreeSE = new TTree(vmsa::vm_tree[0].Data(),vmsa::vm_tree[0].Data());
  mTreeSE = new TTree(SEname.c_str(),SEname.c_str());
  mTreeSE->Branch(vmsa::vm_branch[0].Data(),"StMesonEvent",&mMesonEventSE);
  mTreeSE->SetAutoSave(5000000);

  mMesonEventME = new StMesonEvent();
  std::string MEname = Form("%s_ME",vmsa::vm_tree[0].Data());
  //mTreeME = new TTree(vmsa::vm_tree[0].Data(),vmsa::vm_tree[0].Data());
  mTreeME = new TTree(MEname.c_str(),MEname.c_str());
  mTreeME->Branch(vmsa::vm_branch[0].Data(),"StMesonEvent",&mMesonEventME);
  mTreeME->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------


void StVecMesonTree::WriteMass2Phi()
{
  //h_Mass2->Write();
  //h_KdEdxRig->Write();
  //h_KM2Rig->Write();
  //h_KInvBetaRig->Write(); 
  mTreeSE->Write("",TObject::kOverwrite);
  mTreeME->Write("",TObject::kOverwrite);
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
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

  for(Int_t Bin_Event = 0; Bin_Event < vmsa::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      for(Int_t tag = 0; tag < 2; tag++)
      {
        MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge,Bool_t(tag));
        mHelix[key].clear();
        mMomentum[key].clear();
        mMass2[key].clear();
        mDca[key].clear();
        mNHitsFit[key].clear();
        mNHitsMax[key].clear();
        mDEdx[key].clear();
        mNSigma[key].clear();
      }
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}


//void StVecMesonTree::size_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
//{
//  LOG_INFO << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << Bin_vz << ", Psi2 = " << Bin_Psi2 << endm;
//  LOG_INFO << "Buffer_depth = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endm;
//
//  LOG_INFO << "Size of primaryVertex = " << mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].size() << endm;;
//  LOG_INFO << "Size of refMult       = " << mRefMult[cent9][Bin_vz][Bin_Psi2].size() << endm;;
//  LOG_INFO << "---------------------------------------------------------------------------" << endm;
//
//  for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
//  {
//    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
//    LOG_INFO << "Event Number " << Bin_Event << ":" << endm; 
//    LOG_INFO << "Positive Kaon:" << endm;
//    LOG_INFO << "  Size of Helix        = " << mHelix[key].size() << endm;;
//    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
//    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
//    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
//    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
//    LOG_INFO << "  Size of nSigma       = " << mNSigma[key].size() << endm;
//    LOG_INFO << "  Size of Charge       = " << mCharge[key].size() << endm;
//
//    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
//    LOG_INFO << "Negative Kaon:" << endm;
//    LOG_INFO << "  Size of Helix        = " << mHelix[key].size() << endm;
//    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
//    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
//    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
//    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
//    LOG_INFO << "  Size of nSigma       = " << mNSigma[key].size() << endm;
//    LOG_INFO << "  Size of Charge       = " << mCharge[key].size() << endm;
//    LOG_INFO << "---------------------------------------------------------------------------" << endm;
//  }
//}


void StVecMesonTree::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // event header
      mMesonEventSE->clearTrackList();
      mMesonEventSE->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // QVector
      mMesonEventSE->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // Number of Tracks
      mMesonEventSE->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      mMesonEventSE->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEventSE->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
   
      // start to select phi candidate in a event
      MEKey key_plus_probe   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0,0);
      MEKey key_minus_probe  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1,0);
      MEKey key_plus_tagged  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0,1);
      MEKey key_minus_tagged = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1,1);

      TLorentzVector ltrackA, ltrackB;
      // K+ probe + K- tagged
      for(Int_t n_kminus = 0; n_kminus < mHelix[key_minus_tagged].size(); n_kminus++) // second track loop over K- candidates
      {
        if(mMass2[key_minus_tagged][n_kminus] <= 0.16 || mMass2[key_minus_tagged][n_kminus] >= 0.36) continue;

        StThreeVectorF p_vecB = mHelix[key_minus_tagged][n_kminus].cat(mHelix[key_minus_tagged][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
        p_vecB *= mMomentum[key_minus_tagged][n_kminus];
        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);

        for(Int_t n_kplus = 0; n_kplus < mHelix[key_plus_probe].size(); n_kplus++) // first track loop over K+ candidates
        {
  	  StThreeVectorF p_vecA = mHelix[key_plus_probe][n_kplus].cat(mHelix[key_plus_probe][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
            //cout << "p_vecA = (" << p_vecA.x() << ", " << p_vecA.y() << ", " << p_vecA.z() << ")" << endl;
  	  p_vecA *= mMomentum[key_plus_probe][n_kplus];
  	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);

	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  // fill phi candidate into mTree
	  if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	  {
	    mMesonTrackSE = mMesonEventSE->createTrack();
	    mMesonTrackSE->setMass2A(mMass2[key_plus_probe][n_kplus]); // K+
	    mMesonTrackSE->setMass2B(mMass2[key_minus_tagged][n_kminus]); // K-
	    mMesonTrackSE->setNSigA(mNSigma[key_plus_probe][n_kplus]); // K+
	    mMesonTrackSE->setNSigB(mNSigma[key_minus_tagged][n_kminus]); // K-
	    mMesonTrackSE->setDcaA(mDca[key_plus_probe][n_kplus]); // K+
	    mMesonTrackSE->setDcaB(mDca[key_minus_tagged][n_kminus]); // K-
	    mMesonTrackSE->setNHitsFitA(mNHitsFit[key_plus_probe][n_kplus]); // K+
	    mMesonTrackSE->setNHitsFitB(mNHitsFit[key_minus_tagged][n_kminus]); // K-
	    mMesonTrackSE->setNHitsMaxA(mNHitsMax[key_plus_probe][n_kplus]); // K+
	    mMesonTrackSE->setNHitsMaxB(mNHitsMax[key_minus_tagged][n_kminus]); // K-
	    mMesonTrackSE->setDEdxA(mDEdx[key_plus_probe][n_kplus]); // K+
	    mMesonTrackSE->setDEdxB(mDEdx[key_minus_tagged][n_kminus]); // K-
	    mMesonTrackSE->setTrackA(ltrackA); // K+
	    mMesonTrackSE->setTrackB(ltrackB); // K-
	    mMesonTrackSE->setFlagA(Bin_Event); // K+
	    mMesonTrackSE->setFlagB(Bin_Event); // K-
            mMesonTrackSE->setWhichTagged(1); // 1 means K- is tagged
	    //h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
      // K+ tagged + K- probe
      for(Int_t n_kplus = 0; n_kplus < mHelix[key_plus_tagged].size(); n_kplus++) // first track loop over K+ candidates
      {
        if(mMass2[key_plus_tagged][n_kplus] < 0.16 || mMass2[key_plus_tagged][n_kplus] > 0.36) continue;

	StThreeVectorF p_vecA = mHelix[key_plus_tagged][n_kplus].cat(mHelix[key_plus_tagged][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
        //cout << "p_vecA = (" << p_vecA.x() << ", " << p_vecA.y() << ", " << p_vecA.z() << ")" << endl;
	p_vecA *= mMomentum[key_plus_tagged][n_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);

	for(Int_t n_kminus = 0; n_kminus < mHelix[key_minus_probe].size(); n_kminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix[key_minus_probe][n_kminus].cat(mHelix[key_minus_probe][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_minus_probe][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);
	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  // fill phi candidate into mTree
	  if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	  {
	    mMesonTrackSE = mMesonEventSE->createTrack();
	    mMesonTrackSE->setMass2A(mMass2[key_plus_tagged][n_kplus]); // K+
	    mMesonTrackSE->setMass2B(mMass2[key_minus_probe][n_kminus]); // K-
	    mMesonTrackSE->setNSigA(mNSigma[key_plus_tagged][n_kplus]); // K+
	    mMesonTrackSE->setNSigB(mNSigma[key_minus_probe][n_kminus]); // K-
	    mMesonTrackSE->setDcaA(mDca[key_plus_tagged][n_kplus]); // K+
	    mMesonTrackSE->setDcaB(mDca[key_minus_probe][n_kminus]); // K-
	    mMesonTrackSE->setNHitsFitA(mNHitsFit[key_plus_tagged][n_kplus]); // K+
	    mMesonTrackSE->setNHitsFitB(mNHitsFit[key_minus_probe][n_kminus]); // K-
	    mMesonTrackSE->setNHitsMaxA(mNHitsMax[key_plus_tagged][n_kplus]); // K+
	    mMesonTrackSE->setNHitsMaxB(mNHitsMax[key_minus_probe][n_kminus]); // K-
	    mMesonTrackSE->setDEdxA(mDEdx[key_plus_tagged][n_kplus]); // K+
	    mMesonTrackSE->setDEdxB(mDEdx[key_minus_probe][n_kminus]); // K-
	    mMesonTrackSE->setTrackA(ltrackA); // K+
	    mMesonTrackSE->setTrackB(ltrackB); // K-
	    mMesonTrackSE->setFlagA(Bin_Event); // K+
	    mMesonTrackSE->setFlagB(Bin_Event); // K-
            mMesonTrackSE->setWhichTagged(0); // 0 means K+ is tagged
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTreeSE->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    //cout << "Attempt event mixing" << endl;
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_plus_probe   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0,0);
      MEKey key_A_minus_probe  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1,0);
      MEKey key_A_plus_tagged  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0,1);
      MEKey key_A_minus_tagged = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1,1);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_plus_probe   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0,0);
	MEKey key_B_minus_probe  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1,0);
	MEKey key_B_plus_tagged  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0,1);
	MEKey key_B_minus_tagged = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1,1);

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mMesonEventME->clearTrackList();
	  mMesonEventME->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // QVector
	  mMesonEventME->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // Number of Tracks
	  mMesonEventME->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  mMesonEventME->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEventME->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	TLorentzVector ltrackA, ltrackB;
	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	// Event A K+ tagged + Event B K- probe
	for(Int_t n_kplus = 0; n_kplus < mHelix[key_A_plus_tagged].size(); n_kplus++) // first track loop over K+ candidates from event A
	{
          if(mMass2[key_A_plus_tagged][n_kplus] <= 0.16 || mMass2[key_A_plus_tagged][n_kplus] >= 0.36) continue;
	  StThreeVectorF p_vecA(mHelix[key_A_plus_tagged][n_kplus].cat(mHelix[key_A_plus_tagged][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_plus_tagged][n_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+

	  for(Int_t n_kminus = 0; n_kminus < mHelix[key_B_minus_probe].size(); n_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_minus_probe][n_kminus].cat(mHelix[key_B_minus_probe][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_minus_probe][n_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();
	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mMesonTrackME = mMesonEventME->createTrack();
	      mMesonTrackME->setMass2A(mMass2[key_A_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setMass2B(mMass2[key_B_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setNSigA(mNSigma[key_A_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setNSigB(mNSigma[key_B_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setDcaA(mDca[key_A_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setDcaB(mDca[key_B_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setNHitsFitA(mNHitsFit[key_A_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setNHitsFitB(mNHitsFit[key_B_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setNHitsMaxA(mNHitsMax[key_A_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setNHitsMaxB(mNHitsMax[key_B_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setDEdxA(mDEdx[key_A_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setDEdxB(mDEdx[key_B_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setTrackA(ltrackA); // K+
	      mMesonTrackME->setTrackB(ltrackB); // K-
	      mMesonTrackME->setFlagA(Bin_Event_A); // K+
	      mMesonTrackME->setFlagB(Bin_Event_B); // K-
              mMesonTrackME->setWhichTagged(0); // 0 means K+ is tagged
	      //h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
	// Event A K+ probe + Event B K- tagged
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_B_minus_tagged].size(); n_kminus++) // second track loop over K- candidates from event B
	{
          if(mMass2[key_B_minus_tagged][n_kminus] <= 0.16 || mMass2[key_B_minus_tagged][n_kminus] >= 0.36) continue;
	  StThreeVectorF p_vecB(mHelix[key_B_minus_tagged][n_kminus].cat(mHelix[key_B_minus_tagged][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	  p_vecB *= mMomentum[key_B_minus_tagged][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-
	  for(Int_t n_kplus = 0; n_kplus < mHelix[key_A_plus_probe].size(); n_kplus++) // first track loop over K+ candidates from event A
	  {
	    StThreeVectorF p_vecA(mHelix[key_A_plus_probe][n_kplus].cat(mHelix[key_A_plus_probe][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	    p_vecA *= mMomentum[key_A_plus_probe][n_kplus];
	    ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+


	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();
	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mMesonTrackME = mMesonEventME->createTrack();
	      mMesonTrackME->setMass2A(mMass2[key_A_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setMass2B(mMass2[key_B_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setNSigA(mNSigma[key_A_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setNSigB(mNSigma[key_B_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setDcaA(mDca[key_A_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setDcaB(mDca[key_B_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setNHitsFitA(mNHitsFit[key_A_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setNHitsFitB(mNHitsFit[key_B_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setNHitsMaxA(mNHitsMax[key_A_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setNHitsMaxB(mNHitsMax[key_B_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setDEdxA(mDEdx[key_A_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setDEdxB(mDEdx[key_B_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setTrackA(ltrackA); // K+
	      mMesonTrackME->setTrackB(ltrackB); // K-
	      mMesonTrackME->setFlagA(Bin_Event_A); // K+
	      mMesonTrackME->setFlagB(Bin_Event_B); // K-
              mMesonTrackME->setWhichTagged(1); // 1 means K- is tagged
	      //h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}

	// mix K- candidates from A event with K+ candidates from B event
	// Event A K- tagged + Event B K+ probe
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_A_minus_tagged].size(); n_kminus++) // first track loop over K- candidates from event A
	{
          if(mMass2[key_A_minus_tagged][n_kminus] <= 0.16 || mMass2[key_A_minus_tagged][n_kminus] >= 0.36) continue;
	  StThreeVectorF p_vecA(mHelix[key_A_minus_tagged][n_kminus].cat(mHelix[key_A_minus_tagged][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_minus_tagged][n_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-

	  for(Int_t n_kplus = 0; n_kplus < mHelix[key_B_plus_probe].size(); n_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_plus_probe][n_kplus].cat(mHelix[key_B_plus_probe][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_plus_probe][n_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mMesonTrackME = mMesonEventME->createTrack();
	      mMesonTrackME->setMass2A(mMass2[key_B_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setMass2B(mMass2[key_A_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setNSigA(mNSigma[key_B_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setNSigB(mNSigma[key_A_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setDcaA(mDca[key_B_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setDcaB(mDca[key_A_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setNHitsFitA(mNHitsFit[key_B_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setNHitsFitB(mNHitsFit[key_A_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setNHitsMaxA(mNHitsMax[key_B_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setNHitsMaxB(mNHitsMax[key_A_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setDEdxA(mDEdx[key_B_plus_probe][n_kplus]); // K+
	      mMesonTrackME->setDEdxB(mDEdx[key_A_minus_tagged][n_kminus]); // K-
	      mMesonTrackME->setTrackA(ltrackB); // K+
	      mMesonTrackME->setTrackB(ltrackA); // K-
	      mMesonTrackME->setFlagA(Bin_Event_B); // K+
	      mMesonTrackME->setFlagB(Bin_Event_A); // K-
              mMesonTrackME->setWhichTagged(1); // 1 means K- is tagged
	      //h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
	// Event B K+ tagged + Event A K- probe
	for(Int_t n_kplus = 0; n_kplus < mHelix[key_B_plus_tagged].size(); n_kplus++) // second track loop over K+ candidates from event B
	{
          if(mMass2[key_B_plus_tagged][n_kplus] <= 0.16 || mMass2[key_B_plus_tagged][n_kplus] >= 0.36) continue;
	  StThreeVectorF p_vecB(mHelix[key_B_plus_tagged][n_kplus].cat(mHelix[key_B_plus_tagged][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	  p_vecB *= mMomentum[key_B_plus_tagged][n_kplus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+
	  for(Int_t n_kminus = 0; n_kminus < mHelix[key_A_minus_probe].size(); n_kminus++) // first track loop over K- candidates from event A
	  {
	    StThreeVectorF p_vecA(mHelix[key_A_minus_probe][n_kminus].cat(mHelix[key_A_minus_probe][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	    p_vecA *= mMomentum[key_A_minus_probe][n_kminus];
	    ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-


	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mMesonTrackME = mMesonEventME->createTrack();
	      mMesonTrackME->setMass2A(mMass2[key_B_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setMass2B(mMass2[key_A_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setNSigA(mNSigma[key_B_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setNSigB(mNSigma[key_A_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setDcaA(mDca[key_B_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setDcaB(mDca[key_A_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setNHitsFitA(mNHitsFit[key_B_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setNHitsFitB(mNHitsFit[key_A_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setNHitsMaxA(mNHitsMax[key_B_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setNHitsMaxB(mNHitsMax[key_A_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setDEdxA(mDEdx[key_B_plus_tagged][n_kplus]); // K+
	      mMesonTrackME->setDEdxB(mDEdx[key_A_minus_probe][n_kminus]); // K-
	      mMesonTrackME->setTrackA(ltrackB); // K+
	      mMesonTrackME->setTrackB(ltrackA); // K-
	      mMesonTrackME->setFlagA(Bin_Event_B); // K+
	      mMesonTrackME->setFlagB(Bin_Event_A); // K-
              mMesonTrackME->setWhichTagged(0); // 0 means K+ is tagged
	      //h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
      }
    }
    mTreeME->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------



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
      Float_t pT = track->pMom().Pt();
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
// GAVIN COMMENTING OUT THIS SECTION 07.24.2024 
// The point of commenting this out is to implement a mixed selection of kaons. 
        //cout << "fabs(eta) = " << fabs(eta) << endl;
    	//if( fabs(eta) <= 1.0 && (Mass2 < Mass2_low || Mass2 > Mass2_up) ) continue;
        if( fabs(eta) > 1.5) continue;
// GAVIN COMMENTING OUT THIS SECTION 07.24.2024 
// The point of commenting this out is to implement a mixed selection of kaons. 

    	//if( pT < 0.9 && fabs(eta) <= 1.0 && (Mass2 < Mass2_low || Mass2 > Mass2_up) ) continue;
        //if( fabs(eta) > 1.0) continue;
	
        //cout << "Polarity = " << Polarity << endl;

          

        StThreeVectorD primMom; 
        primMom.setX(track->pMom().x());
        primMom.setY(track->pMom().y());
        primMom.setZ(track->pMom().z());
        StThreeVectorD primVer; 
        primVer.setX(event->primaryVertex().x());
        primVer.setY(event->primaryVertex().y());
        primVer.setZ(event->primaryVertex().z());

        bool tag = gRandom->Uniform(0,1) > 0.5; // 0 if probe (no TOF), 1 if tagged (with m2 TOF cut)
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge,tag);
        
        //mTagged[key].push_back(static_cast<Bool_t>(tag));

	mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
	mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
	mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	mNHitsMax[key].push_back(static_cast<Float_t>(track->nHitsMax())); // nHitsFit
	mDEdx[key].push_back(static_cast<Float_t>(track->dEdx())); // nHitsFit
	mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaKaon()))); // nSigma
	mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
        //cout << "MagneticField = " << event->bField() << "    MAGFIELDFACTOR = " << MAGFIELDFACTOR << "    product = " << event->bField()*MAGFIELDFACTOR << endl; 
	mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
        //h_KdEdxRig->Fill(momentum*Polarity,track->dEdx());
        //h_KM2Rig->Fill(momentum*Polarity,Mass2);
        //h_KInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
	
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  // First we do the SE combinations and this should always be done before any event mixing
  // This is where we will add the TOF tag, which will then be read by the event mixing step:w
  doPhi(0,cent9,Bin_vz,Bin_Psi2);

  if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth)
  {
    doPhi(1,cent9,Bin_vz,Bin_Psi2);
    clear_phi(cent9,Bin_vz,Bin_Psi2);
  }
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
