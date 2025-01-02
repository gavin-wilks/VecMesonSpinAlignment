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
      }
    }
  }
  cout << "Creating a StMesonEvent" << endl;
  mMesonEvent = new StMesonEvent();
  cout << "Creating a TTree" << endl;
  mTree = new TTree(vmsa::vm_tree[0].Data(),vmsa::vm_tree[0].Data());
  cout << "Creating a branch" << endl;
  mTree->Branch(vmsa::vm_branch[0].Data(),"StMesonEvent",&mMesonEvent);
  cout << "SetAutoSave" << endl;
  mTree->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------


void StVecMesonTree::InitRho()
{
  mVecMesonCut = new StVecMesonCut(mEnergy);
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,400,0.5,0.9);
  HistName = "pi_dEdx_Rig";
  h_PidEdxRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,400,0,40);
  HistName = "pi_Mass2_Rig";
  h_PiM2Rig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,100,-0.25,0.5);
  HistName = "pi_InvBeta_Rig";
  h_PiInvBetaRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,300,-0.0,3.0);

  for(Int_t cent = 0; cent < vmsa::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < vmsa::Bin_VertexZ; vz++)
    {
      for(Int_t rho_psi = 0; rho_psi < vmsa::Bin_Rho_Psi; rho_psi++)
      {
        mEventCounter2[cent][vz][rho_psi] = 0;
	clear_rho(cent,vz,rho_psi);
      }
    }
  }

  mMesonEvent = new StMesonEvent();
  mTree = new TTree(vmsa::vm_tree[1].Data(),vmsa::vm_tree[1].Data());
  mTree->Branch(vmsa::vm_branch[1].Data(),"StMesonEvent",&mMesonEvent);
  mTree->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------


void StVecMesonTree::InitKStar()
{
  mVecMesonCut = new StVecMesonCut(mEnergy);
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,400,0.60,1.2);
  HistName = "K_dEdx_Rig";
  h_KdEdxRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,400,0,40);
  HistName = "K_Mass2_Rig";
  h_KM2Rig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,100,-0.25,0.5);
  HistName = "K_InvBeta_Rig";
  h_KInvBetaRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,300,-0.0,3.0);
  HistName = "pi_dEdx_Rig";
  h_PidEdxRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,400,0,40);
  HistName = "pi_Mass2_Rig";
  h_PiM2Rig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,100,-0.25,0.5);
  HistName = "pi_InvBeta_Rig";
  h_PiInvBetaRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,300,-0.0,3.0);

  for(Int_t cent = 0; cent < vmsa::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < vmsa::Bin_VertexZ; vz++)
    {
      for(Int_t kstar_psi = 0; kstar_psi < vmsa::Bin_KStar_Psi; kstar_psi++)
      {
        mEventCounter2[cent][vz][kstar_psi] = 0;
	clear_kstar(cent,vz,kstar_psi);
      }
    }
  }

  mMesonEvent = new StMesonEvent();
  mTree = new TTree(vmsa::vm_tree[2].Data(),vmsa::vm_tree[2].Data());
  mTree->Branch(vmsa::vm_branch[2].Data(),"StMesonEvent",&mMesonEvent);
  mTree->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------



void StVecMesonTree::WriteMass2Phi()
{
  //h_Mass2->Write();
  //h_KdEdxRig->Write();
  //h_KM2Rig->Write();
  //h_KInvBetaRig->Write(); 
  mTree->Write("",TObject::kOverwrite);
}

void StVecMesonTree::WriteMass2Rho()
{
  h_Mass2->Write();
  h_PidEdxRig->Write();
  h_PiM2Rig->Write();
  h_PiInvBetaRig->Write(); 
  mTree->Write("",TObject::kOverwrite);
}

void StVecMesonTree::WriteMass2KStar()
{
  h_Mass2->Write();
  h_KdEdxRig->Write();
  h_KM2Rig->Write();
  h_KInvBetaRig->Write(); 
  h_PidEdxRig->Write();
  h_PiM2Rig->Write();
  h_PiInvBetaRig->Write(); 
  mTree->Write("",TObject::kOverwrite);
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
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
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
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StVecMesonTree::clear_rho(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
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
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
      mHelix[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigma[key].clear();
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StVecMesonTree::clear_kstar(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
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
    for(Int_t id = 0; id < 4; id++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,id);
      mHelix[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigma[key].clear();
      mCharge[key].clear();
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

 void StVecMesonTree::size_rho(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
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
    LOG_INFO << "Positive Pion:" << endm;
    LOG_INFO << "  Size of Helix        = " << mHelix[key].size() << endm;;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigma       = " << mNSigma[key].size() << endm;
    LOG_INFO << "  Size of Charge       = " << mCharge[key].size() << endm;

    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
    LOG_INFO << "Negative Pion:" << endm;
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

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::size_kstar(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
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

    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,2);
    LOG_INFO << "Event Number " << Bin_Event << ":" << endm; 
    LOG_INFO << "Positive Pion:" << endm;
    LOG_INFO << "  Size of Helix        = " << mHelix[key].size() << endm;;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigma       = " << mNSigma[key].size() << endm;
    LOG_INFO << "  Size of Charge       = " << mCharge[key].size() << endm;

    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,3);
    LOG_INFO << "Negative Pion:" << endm;
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

void StVecMesonTree::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // event header
      mMesonEvent->clearTrackList();
      mMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // QVector
      mMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // Number of Tracks
      mMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      mMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
   
      // start to select phi candidate in a event
      MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
      MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);

      TLorentzVector ltrackA, ltrackB;
      for(Int_t n_kplus = 0; n_kplus < mHelix[key_plus].size(); n_kplus++) // first track loop over K+ candidates
      {
	StThreeVectorF p_vecA = mHelix[key_plus][n_kplus].cat(mHelix[key_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
        //cout << "p_vecA = (" << p_vecA.x() << ", " << p_vecA.y() << ", " << p_vecA.z() << ")" << endl;
	p_vecA *= mMomentum[key_plus][n_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);

	for(Int_t n_kminus = 0; n_kminus < mHelix[key_minus].size(); n_kminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix[key_minus][n_kminus].cat(mHelix[key_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_minus][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);
	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  // fill phi candidate into mTree
	  if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	  {
	    mMesonTrack = mMesonEvent->createTrack();
	    mMesonTrack->setMass2A(mMass2[key_plus][n_kplus]); // K+
	    mMesonTrack->setMass2B(mMass2[key_minus][n_kminus]); // K-
	    mMesonTrack->setNSigA(mNSigma[key_plus][n_kplus]); // K+
	    mMesonTrack->setNSigB(mNSigma[key_minus][n_kminus]); // K-
	    mMesonTrack->setDcaA(mDca[key_plus][n_kplus]); // K+
	    mMesonTrack->setDcaB(mDca[key_minus][n_kminus]); // K-
	    mMesonTrack->setNHitsFitA(mNHitsFit[key_plus][n_kplus]); // K+
	    mMesonTrack->setNHitsFitB(mNHitsFit[key_minus][n_kminus]); // K-
	    mMesonTrack->setNHitsMaxA(mNHitsMax[key_plus][n_kplus]); // K+
	    mMesonTrack->setNHitsMaxB(mNHitsMax[key_minus][n_kminus]); // K-
	    mMesonTrack->setDEdxA(mDEdx[key_plus][n_kplus]); // K+
	    mMesonTrack->setDEdxB(mDEdx[key_minus][n_kminus]); // K-
	    mMesonTrack->setTrackA(ltrackA); // K+
	    mMesonTrack->setTrackB(ltrackB); // K-
	    mMesonTrack->setFlagA(Bin_Event); // K+
	    mMesonTrack->setFlagB(Bin_Event); // K-
	    mMesonTrack->setWeightA(mWeight[key_plus][n_kplus]); // K+
	    mMesonTrack->setWeightB(mWeight[key_minus][n_kminus]); // K-
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTree->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    //cout << "Attempt event mixing" << endl;
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
      MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
	MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mMesonEvent->clearTrackList();
	  mMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // QVector
	  mMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // Number of Tracks
	  mMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  mMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	TLorentzVector ltrackA, ltrackB;
	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	for(Int_t n_kplus = 0; n_kplus < mHelix[key_A_plus].size(); n_kplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_plus][n_kplus].cat(mHelix[key_A_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_plus][n_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+

	  for(Int_t n_kminus = 0; n_kminus < mHelix[key_B_minus].size(); n_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_minus][n_kminus].cat(mHelix[key_B_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_minus][n_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();
	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2A(mMass2[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setMass2B(mMass2[key_B_minus][n_kminus]); // K-
	      mMesonTrack->setNSigA(mNSigma[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setNSigB(mNSigma[key_B_minus][n_kminus]); // K-
	      mMesonTrack->setDcaA(mDca[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setDcaB(mDca[key_B_minus][n_kminus]); // K-
	      mMesonTrack->setNHitsFitA(mNHitsFit[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setNHitsFitB(mNHitsFit[key_B_minus][n_kminus]); // K-
	      mMesonTrack->setNHitsMaxA(mNHitsMax[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setNHitsMaxB(mNHitsMax[key_B_minus][n_kminus]); // K-
	      mMesonTrack->setDEdxA(mDEdx[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setDEdxB(mDEdx[key_B_minus][n_kminus]); // K-
	      mMesonTrack->setTrackA(ltrackA); // K+
	      mMesonTrack->setTrackB(ltrackB); // K-
	      mMesonTrack->setFlagA(Bin_Event_A); // K+
	      mMesonTrack->setFlagB(Bin_Event_B); // K-
	      mMesonTrack->setWeightA(mWeight[key_A_plus][n_kplus]); // K+
	      mMesonTrack->setWeightB(mWeight[key_B_minus][n_kminus]); // K-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}

	// mix K- candidates from A event with K+ candidates from B event
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_A_minus].size(); n_kminus++) // first track loop over K- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_minus][n_kminus].cat(mHelix[key_A_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_minus][n_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-

	  for(Int_t n_kplus = 0; n_kplus < mHelix[key_B_plus].size(); n_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_plus][n_kplus].cat(mHelix[key_B_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_plus][n_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2A(mMass2[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setMass2B(mMass2[key_A_minus][n_kminus]); // K-
	      mMesonTrack->setNSigA(mNSigma[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setNSigB(mNSigma[key_A_minus][n_kminus]); // K-
	      mMesonTrack->setDcaA(mDca[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setDcaB(mDca[key_A_minus][n_kminus]); // K-
	      mMesonTrack->setNHitsFitA(mNHitsFit[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setNHitsFitB(mNHitsFit[key_A_minus][n_kminus]); // K-
	      mMesonTrack->setNHitsMaxA(mNHitsMax[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setNHitsMaxB(mNHitsMax[key_A_minus][n_kminus]); // K-
	      mMesonTrack->setDEdxA(mDEdx[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setDEdxB(mDEdx[key_A_minus][n_kminus]); // K-
	      mMesonTrack->setTrackA(ltrackB); // K+
	      mMesonTrack->setTrackB(ltrackA); // K-
	      mMesonTrack->setFlagA(Bin_Event_B); // K+
	      mMesonTrack->setFlagB(Bin_Event_A); // K-
	      mMesonTrack->setWeightA(mWeight[key_B_plus][n_kplus]); // K+
	      mMesonTrack->setWeightB(mWeight[key_A_minus][n_kminus]); // K-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
      }
    }
    mTree->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::doRho(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // event header
      mMesonEvent->clearTrackList();
      mMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // QVector
      mMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // Number of Tracks
      mMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      
      mMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
   
      // start to select phi candidate in a event
      MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
      MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);

      TLorentzVector ltrackA, ltrackB;
      for(Int_t n_piplus = 0; n_piplus < mHelix[key_plus].size(); n_piplus++) // first track loop over K+ candidates
      {
	StThreeVectorF p_vecA = mHelix[key_plus][n_piplus].cat(mHelix[key_plus][n_piplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	p_vecA *= mMomentum[key_plus][n_piplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassPion);
	for(Int_t n_piminus = 0; n_piminus < mHelix[key_minus].size(); n_piminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix[key_minus][n_piminus].cat(mHelix[key_minus][n_piminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_minus][n_piminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassPion);
	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  // fill phi candidate into mTree
	  if(InvMassAB > 0.5 && InvMassAB < 1.1) 
	  {
	    mMesonTrack = mMesonEvent->createTrack();
	    mMesonTrack->setMass2A(mMass2[key_plus][n_piplus]); // K+
	    mMesonTrack->setMass2B(mMass2[key_minus][n_piminus]); // K-
	    mMesonTrack->setNSigA(mNSigma[key_plus][n_piplus]); // K+
	    mMesonTrack->setNSigB(mNSigma[key_minus][n_piminus]); // K-
	    mMesonTrack->setDcaA(mDca[key_plus][n_piplus]); // K+
	    mMesonTrack->setDcaB(mDca[key_minus][n_piminus]); // K-
	    mMesonTrack->setTrackA(ltrackA); // K+
	    mMesonTrack->setTrackB(ltrackB); // K-
	    mMesonTrack->setFlagA(Bin_Event); // K+
	    mMesonTrack->setFlagB(Bin_Event); // K-
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTree->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
      MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
	MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mMesonEvent->clearTrackList();
	  mMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // QVector
	  mMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // Number of Tracks
	  mMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  mMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	TLorentzVector ltrackA, ltrackB;

	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	for(Int_t n_piplus = 0; n_piplus < mHelix[key_A_plus].size(); n_piplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_plus][n_piplus].cat(mHelix[key_A_plus][n_piplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_plus][n_piplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassPion); // K+

	  for(Int_t n_piminus = 0; n_piminus < mHelix[key_B_minus].size(); n_piminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_minus][n_piminus].cat(mHelix[key_B_minus][n_piminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_minus][n_piminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassPion); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > 0.5 && InvMassAB < 1.1) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2A(mMass2[key_A_plus][n_piplus]); // K+
	      mMesonTrack->setMass2B(mMass2[key_B_minus][n_piminus]); // K-
	      mMesonTrack->setNSigA(mNSigma[key_A_plus][n_piplus]); // K+
	      mMesonTrack->setNSigB(mNSigma[key_B_minus][n_piminus]); // K-
	      mMesonTrack->setDcaA(mDca[key_A_plus][n_piplus]); // K+
	      mMesonTrack->setDcaB(mDca[key_B_minus][n_piminus]); // K-
	      mMesonTrack->setTrackA(ltrackA); // K+
	      mMesonTrack->setTrackB(ltrackB); // K-
	      mMesonTrack->setFlagA(Bin_Event_A); // K+
	      mMesonTrack->setFlagB(Bin_Event_B); // K-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}

	// mix K- candidates from A event with K+ candidates from B event
	for(Int_t n_piminus = 0; n_piminus < mHelix[key_A_minus].size(); n_piminus++) // first track loop over K- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_minus][n_piminus].cat(mHelix[key_A_minus][n_piminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_minus][n_piminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassPion); // K-

	  for(Int_t n_piplus = 0; n_piplus < mHelix[key_B_plus].size(); n_piplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_plus][n_piplus].cat(mHelix[key_B_plus][n_piplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_plus][n_piplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassPion); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree
	    if(InvMassAB > 0.5 && InvMassAB < 1.1) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2A(mMass2[key_B_plus][n_piplus]); // K+
	      mMesonTrack->setMass2B(mMass2[key_A_minus][n_piminus]); // K-
	      mMesonTrack->setNSigA(mNSigma[key_B_plus][n_piplus]); // K+
	      mMesonTrack->setNSigB(mNSigma[key_A_minus][n_piminus]); // K-
	      mMesonTrack->setDcaA(mDca[key_B_plus][n_piplus]); // K+
	      mMesonTrack->setDcaB(mDca[key_A_minus][n_piminus]); // K-
	      mMesonTrack->setTrackA(ltrackB); // K+
	      mMesonTrack->setTrackB(ltrackA); // K-
	      mMesonTrack->setFlagA(Bin_Event_B); // K+
	      mMesonTrack->setFlagB(Bin_Event_A); // K-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
      }
    }
    mTree->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------



void StVecMesonTree::doKStar(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // event header
      mMesonEvent->clearTrackList();
      mMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // QVector
      mMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // Number of Tracks
      mMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      mMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
   
      // start to select phi candidate in a event
      MEKey key_kplus   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0); //K+
      MEKey key_kminus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1); //K-
      MEKey key_piplus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,2); //pi+
      MEKey key_piminus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,3); //pi-

      TLorentzVector ltrackA, ltrackB;
      for(Int_t n_kplus = 0; n_kplus < mHelix[key_kplus].size(); n_kplus++) // first track loop over K+ candidates
      {
	StThreeVectorF p_vecA = mHelix[key_kplus][n_kplus].cat(mHelix[key_kplus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	p_vecA *= mMomentum[key_kplus][n_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);
	for(Int_t n_piminus = 0; n_piminus < mHelix[key_piminus].size(); n_piminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix[key_piminus][n_piminus].cat(mHelix[key_piminus][n_piminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_piminus][n_piminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassPion);
	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  if(InvMassAB > 0.60 && InvMassAB < 1.2) 
	  {
	    mMesonTrack = mMesonEvent->createTrack();
	    mMesonTrack->setMass2A(mMass2[key_kplus][n_kplus]); // K+
	    mMesonTrack->setMass2B(mMass2[key_piminus][n_piminus]); // pi-
	    mMesonTrack->setNSigA(mNSigma[key_kplus][n_kplus]); // K+
	    mMesonTrack->setNSigB(mNSigma[key_piminus][n_piminus]); // pi-
	    mMesonTrack->setDcaA(mDca[key_kplus][n_kplus]); // K+
	    mMesonTrack->setDcaB(mDca[key_piminus][n_piminus]); // pi-
	    mMesonTrack->setChargeA(mCharge[key_kplus][n_kplus]); // K+
	    mMesonTrack->setChargeB(mCharge[key_piminus][n_piminus]); // pi-
	    mMesonTrack->setTrackA(ltrackA); // K+
	    mMesonTrack->setTrackB(ltrackB); // pi-
	    mMesonTrack->setFlagA(Bin_Event); // K+
	    mMesonTrack->setFlagB(Bin_Event); // pi-
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }

      for(Int_t n_piplus = 0; n_piplus < mHelix[key_piplus].size(); n_piplus++) // first track loop over pi+ candidates
      {
	StThreeVectorF p_vecA = mHelix[key_piplus][n_piplus].cat(mHelix[key_piplus][n_piplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	p_vecA *= mMomentum[key_piplus][n_piplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassPion);
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_kminus].size(); n_kminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix[key_kminus][n_kminus].cat(mHelix[key_kminus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_kminus][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);
	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  if(InvMassAB > 0.6 && InvMassAB < 1.2) 
	  {
	    mMesonTrack = mMesonEvent->createTrack();
	    mMesonTrack->setMass2B(mMass2[key_piplus][n_piplus]); // pi+
	    mMesonTrack->setMass2A(mMass2[key_kminus][n_kminus]); // K-
	    mMesonTrack->setNSigB(mNSigma[key_piplus][n_piplus]); // pi+
	    mMesonTrack->setNSigA(mNSigma[key_kminus][n_kminus]); // K-
	    mMesonTrack->setDcaB(mDca[key_piplus][n_piplus]); // pi+
	    mMesonTrack->setDcaA(mDca[key_kminus][n_kminus]); // K-
	    mMesonTrack->setChargeB(mCharge[key_piplus][n_piplus]); // pi+
	    mMesonTrack->setChargeA(mCharge[key_kminus][n_kminus]); // K-
	    mMesonTrack->setTrackB(ltrackA); // pi+
	    mMesonTrack->setTrackA(ltrackB); // K-
	    mMesonTrack->setFlagB(Bin_Event); // pi+
	    mMesonTrack->setFlagA(Bin_Event); // K-
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTree->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_kplus   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0); //K+
      MEKey key_A_kminus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1); //K-
      MEKey key_A_piplus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,2); //pi+
      MEKey key_A_piminus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,3); //pi-
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
        MEKey key_B_kplus   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0); //K+
        MEKey key_B_kminus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1); //K-
        MEKey key_B_piplus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,2); //pi+
        MEKey key_B_piminus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,3); //pi-

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mMesonEvent->clearTrackList();
	  mMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // QVector
	  mMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // Number of Tracks
	  mMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  mMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	TLorentzVector ltrackA, ltrackB;

	// start to mix events
	// mix K+ candidates from A event with pi- candidates from B event
	for(Int_t n_kplus = 0; n_kplus < mHelix[key_A_kplus].size(); n_kplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_kplus][n_kplus].cat(mHelix[key_A_kplus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_kplus][n_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+

	  for(Int_t n_piminus = 0; n_piminus < mHelix[key_B_piminus].size(); n_piminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_piminus][n_piminus].cat(mHelix[key_B_piminus][n_piminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_piminus][n_piminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassPion); // pi-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > 0.60 && InvMassAB < 1.2) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2A(mMass2[key_A_kplus][n_kplus]); // K+
	      mMesonTrack->setMass2B(mMass2[key_B_piminus][n_piminus]); // pi-
	      mMesonTrack->setNSigA(mNSigma[key_A_kplus][n_kplus]); // K+
	      mMesonTrack->setNSigB(mNSigma[key_B_piminus][n_piminus]); // pi-
	      mMesonTrack->setDcaA(mDca[key_A_kplus][n_kplus]); // K+
	      mMesonTrack->setDcaB(mDca[key_B_piminus][n_piminus]); // pi-
	      mMesonTrack->setChargeA(mCharge[key_A_kplus][n_kplus]); // K+
	      mMesonTrack->setChargeB(mCharge[key_B_piminus][n_piminus]); // pi-
	      mMesonTrack->setTrackA(ltrackA); // K+
	      mMesonTrack->setTrackB(ltrackB); // pi-
	      mMesonTrack->setFlagA(Bin_Event_A); // K+
	      mMesonTrack->setFlagB(Bin_Event_B); // pi-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}

        // mix pi+ candidates from A event with K- candidates from B event
	for(Int_t n_piplus = 0; n_piplus < mHelix[key_A_piplus].size(); n_piplus++) // first track loop over pi+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_piplus][n_piplus].cat(mHelix[key_A_piplus][n_piplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_piplus][n_piplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassPion); // pi+

	  for(Int_t n_kminus = 0; n_kminus < mHelix[key_B_kminus].size(); n_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_kminus][n_kminus].cat(mHelix[key_B_kminus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_kminus][n_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > 0.6 && InvMassAB < 1.2) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2B(mMass2[key_A_piplus][n_piplus]); // pi+
	      mMesonTrack->setMass2A(mMass2[key_B_kminus][n_kminus]); // K-
	      mMesonTrack->setNSigB(mNSigma[key_A_piplus][n_piplus]); // pi+
	      mMesonTrack->setNSigA(mNSigma[key_B_kminus][n_kminus]); // K-
	      mMesonTrack->setDcaB(mDca[key_A_piplus][n_piplus]); // pi+
	      mMesonTrack->setDcaA(mDca[key_B_kminus][n_kminus]); // K-
	      mMesonTrack->setChargeB(mCharge[key_A_piplus][n_piplus]); // pi+
	      mMesonTrack->setChargeA(mCharge[key_B_kminus][n_kminus]); // K-
	      mMesonTrack->setTrackB(ltrackA); // pi+
	      mMesonTrack->setTrackA(ltrackB); // K-
	      mMesonTrack->setFlagB(Bin_Event_A); // pi+
	      mMesonTrack->setFlagA(Bin_Event_B); // K-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
	// mix K- candidates from A event with pi+ candidates from B event
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_A_kminus].size(); n_kminus++) // first track loop over K- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_kminus][n_kminus].cat(mHelix[key_A_kminus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_kminus][n_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-

	  for(Int_t n_piplus = 0; n_piplus < mHelix[key_B_piplus].size(); n_piplus++) // second track loop over pi+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_piplus][n_piplus].cat(mHelix[key_B_piplus][n_piplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_piplus][n_piplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassPion); // pi+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > 0.60 && InvMassAB < 1.2) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2B(mMass2[key_B_piplus][n_piplus]); // pi+
	      mMesonTrack->setMass2A(mMass2[key_A_kminus][n_kminus]); // K-
	      mMesonTrack->setNSigB(mNSigma[key_B_piplus][n_piplus]); // pi+
	      mMesonTrack->setNSigA(mNSigma[key_A_kminus][n_kminus]); // K-
	      mMesonTrack->setDcaB(mDca[key_B_piplus][n_piplus]); // pi+
	      mMesonTrack->setDcaA(mDca[key_A_kminus][n_kminus]); // K-
	      mMesonTrack->setChargeB(mCharge[key_B_piplus][n_piplus]); // pi+
	      mMesonTrack->setChargeA(mCharge[key_A_kminus][n_kminus]); // K-
	      mMesonTrack->setTrackB(ltrackB); // pi+
	      mMesonTrack->setTrackA(ltrackA); // K-
	      mMesonTrack->setFlagB(Bin_Event_B); // pi+
	      mMesonTrack->setFlagA(Bin_Event_A); // K-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
        // mix pi- candidates from A event with K+ candidates from B event
	for(Int_t n_piminus = 0; n_piminus < mHelix[key_A_piminus].size(); n_piminus++) // first track loop over pi- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_piminus][n_piminus].cat(mHelix[key_A_piminus][n_piminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_piminus][n_piminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassPion); // pi-

	  for(Int_t n_kplus = 0; n_kplus < mHelix[key_B_kplus].size(); n_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_kplus][n_kplus].cat(mHelix[key_B_kplus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_kplus][n_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > 0.6 && InvMassAB < 1.2) 
	    {
	      mMesonTrack = mMesonEvent->createTrack();
	      mMesonTrack->setMass2A(mMass2[key_B_kplus][n_kplus]); // K+
	      mMesonTrack->setMass2B(mMass2[key_A_piminus][n_piminus]); // pi-
	      mMesonTrack->setNSigA(mNSigma[key_B_kplus][n_kplus]); // K+
	      mMesonTrack->setNSigB(mNSigma[key_A_piminus][n_piminus]); // pi-
	      mMesonTrack->setDcaA(mDca[key_B_kplus][n_kplus]); // K+
	      mMesonTrack->setDcaB(mDca[key_A_piminus][n_piminus]); // pi-
              mMesonTrack->setChargeA(mCharge[key_B_kplus][n_kplus]); // K+
	      mMesonTrack->setChargeB(mCharge[key_A_piminus][n_piminus]); // pi-
	      mMesonTrack->setTrackA(ltrackB); // K+
	      mMesonTrack->setTrackB(ltrackA); // pi-
	      mMesonTrack->setFlagA(Bin_Event_B); // K+
	      mMesonTrack->setFlagB(Bin_Event_A); // pi-
	      h_Mass2->Fill(pt,InvMassAB);
	    }
	  }
	}
      }
    }
    mTree->Fill();
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
  mReweight[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(mReweightFactor));

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
    	if( fabs(eta) <= 1.0 && (Mass2 < Mass2_low || Mass2 > Mass2_up) ) continue;
        if( fabs(eta) >  1.0) continue;
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
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
	mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
	mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	mNHitsMax[key].push_back(static_cast<Float_t>(track->nHitsMax())); // nHitsFit
	mDEdx[key].push_back(static_cast<Float_t>(track->dEdx())); // nHitsFit
	mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaKaon()))); // nSigma
	mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
        mWeight[key].push_back(static_cast<Float_t>(mReweight[cent9][Bin_vz][Bin_Psi2][Bin_Event]));
        //cout << "MagneticField = " << event->bField() << "    MAGFIELDFACTOR = " << MAGFIELDFACTOR << "    product = " << event->bField()*MAGFIELDFACTOR << endl; 
	mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
        //h_KdEdxRig->Fill(momentum*Polarity,track->dEdx());
        //h_KM2Rig->Fill(momentum*Polarity,Mass2);
        //h_KInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
	
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_phi(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth)
    {
      doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_phi(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::MixEvent_Rho(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2)
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz, Bin_Psi2;

  Float_t vz_start = vmsa::mVzMaxMap[mEnergy];
  Float_t vz_bin = 2*vz_start/vmsa::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/vmsa::Bin_Rho_Psi;

  for(Int_t i = 0; i < vmsa::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < vmsa::Bin_Rho_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Enent Information
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

    if(mVecMesonCut->passTrackMeson(track,event,1))
    {
      Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track,pico);
      //Float_t scale_nSigma_factor = vmsa::mSigScaleMap[mEnergy];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().Mag();
      Float_t Mass2_low = -0.20;
      Float_t Mass2_up = 0.15;

      Int_t charge = 0; // pi+
      if(Polarity < 0) charge = 1; // pi-

      if(mVecMesonCut->passSigPionCut(track,1))
      {
	if(
	    (momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
	    || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)
	  )
	{
          StThreeVectorD primMom; 
          primMom.setX(track->pMom().x());
          primMom.setY(track->pMom().y());
          primMom.setZ(track->pMom().z());
          StThreeVectorD primVer; 
          primVer.setX(event->primaryVertex().x());
          primVer.setY(event->primaryVertex().y());
          primVer.setZ(event->primaryVertex().z());
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaPion()))); // nSigma
	  mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
	  mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
          //cout << "Fill track information" << endl;
          //cout << "dEdx = " << track->dEdx() << "   Mass2 = " << Mass2 << "   InvBeta = " << 1.0/mVecMesonCut->getBeta(track,pico) << endl; 
          h_PidEdxRig->Fill(momentum*Polarity,track->dEdx());
          h_PiM2Rig->Fill(momentum*Polarity,Mass2);
          h_PiInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
	}
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doRho(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_rho(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth)
    {
      doRho(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_rho(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------


void StVecMesonTree::MixEvent_KStar(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2)
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz, Bin_Psi2;

  Float_t vz_start = vmsa::mVzMaxMap[mEnergy];
  Float_t vz_bin = 2*vz_start/vmsa::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/vmsa::Bin_KStar_Psi;

  for(Int_t i = 0; i < vmsa::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < vmsa::Bin_KStar_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Enent Information
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

    if(mVecMesonCut->passTrackMeson(track,event,2))
    {
      Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track,pico);
      //Float_t scale_nSigma_factor = vmsa::mSigScaleMap[mEnergy];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().Mag();
      Float_t Mass2K_low = 0.16;
      Float_t Mass2K_up = 0.36;
      Float_t Mass2Pi_low = -0.20;
      Float_t Mass2Pi_up = 0.15;
      
      Int_t charge = 0; // +
      if(Polarity < 0) charge = 1; // -
     
      {
	if(
            (Mass2 > Mass2K_low && Mass2 < Mass2K_up) || (Mass2 < -10.0 && mVecMesonCut->passSigKaonCut(track,2))
	  )
	{
          StThreeVectorD primMom; 
          primMom.setX(track->pMom().x());
          primMom.setY(track->pMom().y());
          primMom.setZ(track->pMom().z());
          StThreeVectorD primVer; 
          primVer.setX(event->primaryVertex().x());
          primVer.setY(event->primaryVertex().y());
          primVer.setZ(event->primaryVertex().z());
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaKaon()))); // nSigma
	  mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
	  mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
          cout << "MagneticField = " << event->bField() << "    MAGFIELDFACTOR = " << MAGFIELDFACTOR << endl; 
          //cout << "Fill track information" << endl;
          //cout << "dEdx = " << track->dEdx() << "   Mass2 = " << Mass2 << "   InvBeta = " << 1.0/mVecMesonCut->getBeta(track,pico) << endl; 
          mCharge[key].push_back(static_cast<Int_t>(charge));
          h_KdEdxRig->Fill(momentum*Polarity,track->dEdx());
          h_KM2Rig->Fill(momentum*Polarity,Mass2);
          h_KInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
	}
      }
      {
        charge += 2; // changes id to pi+ and pi-
	if(
            (Mass2 > Mass2Pi_low && Mass2 < Mass2Pi_up) || (Mass2 < -10.0 && mVecMesonCut->passSigPionCut(track,2)) 
	  )
	{
          StThreeVectorD primMom; 
          primMom.setX(track->pMom().x());
          primMom.setY(track->pMom().y());
          primMom.setZ(track->pMom().z());
          StThreeVectorD primVer; 
          primVer.setX(event->primaryVertex().x());
          primVer.setY(event->primaryVertex().y());
          primVer.setZ(event->primaryVertex().z());
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaPion()))); // nSigma
	  mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
	  mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
          //cout << "Fill track information" << endl;
          //cout << "dEdx = " << track->dEdx() << "   Mass2 = " << Mass2 << "   InvBeta = " << 1.0/mVecMesonCut->getBeta(track,pico) << endl; 
          mCharge[key].push_back(static_cast<Int_t>(charge));
          
          h_PidEdxRig->Fill(momentum*Polarity,track->dEdx());
          h_PiM2Rig->Fill(momentum*Polarity,Mass2);
          h_PiInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
	}
      }

    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doKStar(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_kstar(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth_KStar)
    {
      doKStar(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_kstar(cent9,Bin_vz,Bin_Psi2);
    }
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

void StVecMesonTree::passEvent(Int_t N_prim, Int_t N_non_prim, Int_t N_Tof_match, Float_t reweight)
{
  mNumber_prim = N_prim;
  mNumber_non_prim = N_non_prim;
  mNumber_Tof_match = N_Tof_match;
  mReweightFactor = reweight;
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
