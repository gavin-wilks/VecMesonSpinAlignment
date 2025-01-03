#include "StRoot/StVecMesonMaker/StVecMesonTree.h"
#include "StRoot/StVecMesonMaker/StUtility.h"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/StVecMesonAna/StVecMesonHist.h"
#include "StRoot/StVecMesonAna/StVecMesonHistoFlow.h"
#include "StRoot/StVecMesonAna/StVecMesonFillCorr.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StMesonEvent/StMesonEvent.h"
#include <vector>
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StRoot/StVecMesonMaker/StVecMesonMEKey.h"
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
StVecMesonTree::StVecMesonTree(Int_t energy, Int_t flag)
{
  mEnergy = energy;
  mX_flag = flag;
}

StVecMesonTree::~StVecMesonTree()
{
  /* */
}

//------------------------------------------------------------------------------------------------------------------


/*void StVecMesonTree::InitRho()
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
  mTree = new TTree("RhoMesonEvent","RhoMesonEvent");
  mTree->Branch("rho_SpinAlignment_branch","StMesonEvent",&mMesonEvent);
  //mTree->SetAutoSave(5000000);
}
*/
//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::InitPhi()
{
  mVecMesonCorr = new StVecMesonFillCorr(mEnergy);
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonHistoManger = new StVecMesonHist();  
  mVecMesonHistoFlow = new StVecMesonHistoFlow();  

  mUtility = new StUtility(mEnergy);
  //mRefMultCorr = new StRefMultCorr("refmult");

  mVecMesonCorr->InitReCenterCorrection();
  mVecMesonCorr->InitShiftCorrection();
  mVecMesonCorr->InitResolutionCorr();
  mVecMesonHistoManger->InitSys(mX_flag,0);
  //mVecMesonHistoFlow->InitSys(mX_flag,0);

  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,400,0.60,1.2);
  HistName = "K_dEdx_Rig";
  h_KdEdxRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,400,0,40);
  //HistName = "K_Mass2_Rig";
  //h_KM2Rig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,100,-0.25,0.5);
  //HistName = "K_InvBeta_Rig";
  //h_KInvBetaRig = new TH2F(HistName.Data(),HistName.Data(),400,-4.5,4.5,300,-0.0,3.0);

  for(Int_t cent = 0; cent < vmsa::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < vmsa::Bin_VertexZ; vz++)
    {
      for(Int_t psi = 0; psi < vmsa::Bin_Phi_Psi; psi++)
      {
        mEventCounter2[cent][vz][psi] = 0;
	clear_phi(cent,vz,psi,mX_flag);
      }
    }
  }  
}

void StVecMesonTree::InitKStar()
{
  mVecMesonCorr = new StVecMesonFillCorr(mEnergy);
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonHistoManger = new StVecMesonHist();  
  mVecMesonHistoFlow = new StVecMesonHistoFlow();  

  mUtility = new StUtility(mEnergy);
  //mRefMultCorr = new StRefMultCorr("refmult");

  mVecMesonCorr->InitReCenterCorrection();
  mVecMesonCorr->InitShiftCorrection();
  mVecMesonCorr->InitResolutionCorr();
  mVecMesonHistoManger->InitSys(mX_flag,2);
  mVecMesonHistoFlow->InitSys(mX_flag,2);

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
	clear_kstar(cent,vz,kstar_psi,mX_flag);
      }
    }
  }  
}

//------------------------------------------------------------------------------------------------------------------



/*void StVecMesonTree::WriteMass2Phi()
{
  h_Mass2->Write();
  h_KdEdxRig->Write();
  h_KM2Rig->Write();
  h_KInvBetaRig->Write(); 
  mTree->Write("",TObject::kOverwrite);
}

void StVecMesonTree::WriteMass2Rho()
{
  h_Mass2->Write();
  h_PidEdxRig->Write();
  h_PiM2Rig->Write();
  h_PiInvBetaRig->Write(); 
  mTree->Write("",TObject::kOverwrite);
}*/

void StVecMesonTree::WriteMass2Phi()
{
  h_Mass2->Write();
  h_KdEdxRig->Write();
  //h_KM2Rig->Write();
  //h_KInvBetaRig->Write(); 
  mVecMesonHistoManger->WriteSys(mX_flag,0); 
  //mVecMesonHistoFlow->WriteSys(mX_flag,0); 
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
  mVecMesonHistoManger->WriteSys(mX_flag,2); 
  mVecMesonHistoFlow->WriteSys(mX_flag,2); 
  //mTree->Write("",TObject::kOverwrite);
}
//------------------------------------------------------------------------------------------------------------------

/*void StVecMesonTree::clear_rho(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
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
}*/

void StVecMesonTree::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2, Int_t SM_flag)
{
  int numEv = 0;
  if(SM_flag == 0) numEv = 1;
  if(SM_flag == 0) numEv = vmsa::Buffer_depth;

  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  //mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  //mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunIdx[cent9][Bin_vz][Bin_Psi2].clear();
  mReweight[cent9][Bin_vz][Bin_Psi2].clear();
  mVzEPBin[cent9][Bin_vz][Bin_Psi2].clear();
  //mEventId[cent9][Bin_vz][Bin_Psi2].clear();
  //mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  //mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  //mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  //mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  //mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2East[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2West[cent9][Bin_vz][Bin_Psi2].clear();
  //mQ2Full[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].clear();
  //mNumTrackFull[cent9][Bin_vz][Bin_Psi2].clear();
  //mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].clear();
  //mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t Bin_Event = 0; Bin_Event < numEv; Bin_Event++)
  {
    for(Int_t id = 0; id < 2; id++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,id);
      mHelix[key].clear();
      mMomentum[key].clear();
      //mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigma[key].clear();
      mCharge[key].clear();
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}


void StVecMesonTree::clear_kstar(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2, Int_t SM_flag)
{
  int numEv = 0;
  if(SM_flag == 0) numEv = 1;
  if(SM_flag == 0) numEv = vmsa::Buffer_depth;

  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  //mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  //mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunIdx[cent9][Bin_vz][Bin_Psi2].clear();
  mReweight[cent9][Bin_vz][Bin_Psi2].clear();
  mVzEPBin[cent9][Bin_vz][Bin_Psi2].clear();
  //mEventId[cent9][Bin_vz][Bin_Psi2].clear();
  //mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  //mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  //mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  //mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  //mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2East[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2West[cent9][Bin_vz][Bin_Psi2].clear();
  //mQ2Full[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].clear();
  //mNumTrackFull[cent9][Bin_vz][Bin_Psi2].clear();
  //mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].clear();
  //mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t Bin_Event = 0; Bin_Event < numEv; Bin_Event++)
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

/*
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
*/
//------------------------------------------------------------------------------------------------------------------

/*void StVecMesonTree::size_kstar(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
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
}*/

//------------------------------------------------------------------------------------------------------------------

/*void StVecMesonTree::doRho(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
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
	  if(InvMassAB > vmsa::mMassPion*2 && InvMassAB < 0.9) 
	  {
	    mesonTrack = mMesonEvent->createTrack();
	    mesonTrack->setMass2A(mMass2[key_plus][n_piplus]); // K+
	    mesonTrack->setMass2B(mMass2[key_minus][n_piminus]); // K-
	    mesonTrack->setNSigA(mNSigma[key_plus][n_piplus]); // K+
	    mesonTrack->setNSigB(mNSigma[key_minus][n_piminus]); // K-
	    mesonTrack->setDcaA(mDca[key_plus][n_piplus]); // K+
	    mesonTrack->setDcaB(mDca[key_minus][n_piminus]); // K-
	    mesonTrack->setTrackA(ltrackA); // K+
	    mesonTrack->setTrackB(ltrackB); // K-
	    mesonTrack->setFlagA(Bin_Event); // K+
	    mesonTrack->setFlagB(Bin_Event); // K-
	  }

	  // Fill histogram with InvMassAB information
	  h_Mass2->Fill(pt,InvMassAB);
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

	    // fill phi candidate background into mTree
	    if(InvMassAB > vmsa::mMassPion*2 && InvMassAB < 1.05) 
	    {
	      mesonTrack = mMesonEvent->createTrack();
	      mesonTrack->setMass2A(mMass2[key_A_plus][n_piplus]); // K+
	      mesonTrack->setMass2B(mMass2[key_B_minus][n_piminus]); // K-
	      mesonTrack->setNSigA(mNSigma[key_A_plus][n_piplus]); // K+
	      mesonTrack->setNSigB(mNSigma[key_B_minus][n_piminus]); // K-
	      mesonTrack->setDcaA(mDca[key_A_plus][n_piplus]); // K+
	      mesonTrack->setDcaB(mDca[key_B_minus][n_piminus]); // K-
	      mesonTrack->setTrackA(ltrackA); // K+
	      mesonTrack->setTrackB(ltrackB); // K-
	      mesonTrack->setFlagA(Bin_Event_A); // K+
	      mesonTrack->setFlagB(Bin_Event_B); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
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
	    if(InvMassAB > vmsa::mMassPion*2 && InvMassAB < 1.05) 
	    {
	      mesonTrack = mMesonEvent->createTrack();
	      mesonTrack->setMass2A(mMass2[key_B_plus][n_piplus]); // K+
	      mesonTrack->setMass2B(mMass2[key_A_minus][n_piminus]); // K-
	      mesonTrack->setNSigA(mNSigma[key_B_plus][n_piplus]); // K+
	      mesonTrack->setNSigB(mNSigma[key_A_minus][n_piminus]); // K-
	      mesonTrack->setDcaA(mDca[key_B_plus][n_piplus]); // K+
	      mesonTrack->setDcaB(mDca[key_A_minus][n_piminus]); // K-
	      mesonTrack->setTrackA(ltrackB); // K+
	      mesonTrack->setTrackB(ltrackA); // K-
	      mesonTrack->setFlagA(Bin_Event_B); // K+
	      mesonTrack->setFlagB(Bin_Event_A); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTree->Fill();
  }
}
*/
//------------------------------------------------------------------------------------------------------------------


void StVecMesonTree::FillPhi(StMesonTrack* Meson_track, Int_t mCent, Int_t vz, Int_t Psi2)
{

  // Initialise Track 
  Float_t m2A = -999.9;
  Float_t m2B = -999.9;
  Float_t nsA = -999.9;
  Float_t nsB = -999.9;
  Float_t dcaA = -999.9;
  Float_t dcaB = -999.9;
  Float_t nhA = -1;
  Float_t nhB = -1;
  TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
  TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
  Int_t flagA = -1;
  Int_t flagB = -1;

  TVector2 mQv2West = mQ2West[mCent][vz][Psi2][0];
  TVector2 mQv2East = mQ2East[mCent][vz][Psi2][0];

  int mVzSign = mVzEPBin[mCent][vz][Psi2][0];
  int mRunIndex = mRunIdx[mCent][vz][Psi2][0];
  int reweight = mReweight[mCent][vz][Psi2][0];

  // get Track Information
  if(mVecMesonCorr->passTrackEtaNumCut(mNumTrackEast[mCent][vz][Psi2][0],mNumTrackWest[mCent][vz][Psi2][0]))
  {
    //m2A = Meson_track->getMass2A();
    //m2B = Meson_track->getMass2B();
    nsA = Meson_track->getNSigA();
    nsB = Meson_track->getNSigB();
    dcaA = Meson_track->getDcaA();
    dcaB = Meson_track->getDcaB();
    nhA = Meson_track->getNHitA();
    nhB = Meson_track->getNHitB();
    lTrackA = Meson_track->getTrackA();
    lTrackB = Meson_track->getTrackB();
    flagA = Meson_track->getFlagA();
    flagB = Meson_track->getFlagB();
    Float_t pA = lTrackA.P();
    Float_t pB = lTrackB.P();
    TLorentzVector lTrack = lTrackA + lTrackB; // Phi meson
    Float_t pt_lTrack = lTrack.Perp();
    Float_t rapidity_lTrack = lTrack.Rapidity();
    //cout << "Rapidity = " << TMath::Abs(rapidity_lTrack) << endl;
    if(TMath::Abs(rapidity_lTrack) > vmsa::mRapidityMax) return;
    //cout << "Pass rapidity cut" << endl;

    Float_t InvMass_lTrack = lTrack.M();
    TVector3 vBeta = -1.0*lTrack.BoostVector(); // get phi beta
    TLorentzVector lKRest = lTrackA;
    lKRest.Boost(vBeta); // boost kaon back to phi rest frame
    TVector3 vKRest = lKRest.Vect().Unit(); // kaon momentum direction in phi rest frame
    
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++) // systematic loop for dca
    {
      if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,0)) ) continue;
      mVecMesonHistoManger->FillDcaSys(dcaA,dcaB,i_dca); // fill QA for dcaA and dcaB

      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++) // systematic loop for nSigma
      {
        if(i_dca != 0 && i_sig != 0) continue;
        if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,0)) ) continue;
        mVecMesonHistoManger->FillSigSys(nsA,nsB,i_sig); // fill QA for nsA and nsB 
    
        //for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++) // systematic loop for nSigma
        //{
          //if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
          //if( !(mVecMesonCut->passTrackNHitSys(nhA,nhB,i_nhit,0)) ) continue;
          //mVecMesonHistoManger->FillNHitSys(nhA,nhB,i_nhit); // fill QA for nhA and nhB 

          //////////////// GLOBAL SPIN ALIGNMENT ///////////////////////////////////////////////////
          if(mVecMesonCut->passEtaEast(lTrackA)) // K+/- neg eta(east)
          { // Below is West Only
            TVector2 Q2Vector = mQv2West;
            // subtract auto-correlation from pos eta(west) event plane
            if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackB);
              TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
              TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
            }
            Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
            Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,mRunIndex,mCent,mVzSign);

            TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
            TVector3 nQ = nQ_West.Unit();
            Double_t CosThetaStar = vKRest.Dot(nQ);
       
            //cout << "Fill Sys" << endl;
            mVecMesonHistoManger->FillSys(pt_lTrack,mCent,CosThetaStar,i_dca,i_sig,0,Res2,InvMass_lTrack,reweight,mX_flag,0);
   
            //TVector3 nQ_West_EP(TMath::Cos(Psi2_west),TMath::Sin(Psi2_west),0.0); // tangent vector of 2nd Event Plane
            //TVector3 nQ_EP = nQ_West_EP.Unit();          
            //Double_t CosThetaStar_EP = vKRest.Dot(nQ_EP);
          }
    
          if(mVecMesonCut->passEtaWest(lTrackA)) // K+/- pos eta (west)
          { // Below is East Only
            TVector2 Q2Vector = mQv2East;
            // subtract auto-correlation from pos eta(west) event plane
            if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackB);
              TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
              TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
            }
            Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
            Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,mRunIndex,mCent,mVzSign);
            
            TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane
            TVector3 nQ = nQ_East.Unit();
            Double_t CosThetaStar = vKRest.Dot(nQ);

            //cout << "Fill Sys" << endl;
            mVecMesonHistoManger->FillSys(pt_lTrack,mCent,CosThetaStar,i_dca,i_sig,0,Res2,InvMass_lTrack,reweight,mX_flag,0);
            
            //TVector3 nQ_East_EP(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0); // tangent vector of 2nd Event Plane
            //TVector3 nQ_EP = nQ_East_EP.Unit();
            //Double_t CosThetaStar_EP = vKRest.Dot(nQ_EP);
            //mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,2);
          }
          //////////////// GLOBAL SPIN ALIGNMENT ///////////////////////////////////////////////////
          
          //////////////// ELLIPTIC FLOW ///////////////////////////////////////////////////
          //if(mVecMesonCut->passEtaEast(lTrack)) // phi neg eta(east)
          //{ // Below is West Only
          //  TVector2 Q2Vector = mQv2West;
          //  // subtract auto-correlation from pos eta(west) event plane
          //  if(flagA == 0 && mVecMesonCut->passTrackEP(lTrackA,dcaA) && mVecMesonCut->passTrackEtaWest(lTrackA)) // trackB
          //  {
          //    Float_t  w = mVecMesonCorr->getWeight(lTrackA);
          //    TVector2 q2VectorA = mVecMesonCorr->calq2Vector(lTrackA);
          //    TVector2 q2CorrA   = mVecMesonCorr->getReCenterPar_West(mCent,mRunIndex,mVzSign);
          //    Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);
          //  }
          //  if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
          //  {
          //    Float_t  w = mVecMesonCorr->getWeight(lTrackB);
          //    TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
          //    TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(mCent,mRunIndex,mVzSign);
          //    Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
          //  }
          //  Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
          //  Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,mRunIndex,mCent,mVzSign);
          //  Float_t PhiMinusPsi = lTrack.Phi() - Psi2_west;
          //  
          //  if(PhiMinusPsi >= -(3.0/2.0)*TMath::Pi() && PhiMinusPsi < -(1.0/2.0)*TMath::Pi()) PhiMinusPsi += TMath::Pi();
          //  if(PhiMinusPsi <=  (3.0/2.0)*TMath::Pi() && PhiMinusPsi >  (1.0/2.0)*TMath::Pi()) PhiMinusPsi -= TMath::Pi();
          //  PhiMinusPsi = TMath::Abs(PhiMinusPsi);
          // 
          //  mVecMesonHistoFlow->FillSys(pt_lTrack,mCent,PhiMinusPsi,i_dca,i_sig,i_nhit,Res2,InvMass_lTrack,reweight,mX_flag,2);
          //}
  
          //if(mVecMesonCut->passEtaWest(lTrack)) // phi pos eta (west)
          //{ // Below is East Only
          //  TVector2 Q2Vector = mQv2East;
          //  // subtract auto-correlation from pos eta(west) event plane
          //  if(flagA == 0 && mVecMesonCut->passTrackEP(lTrackA,dcaA) && mVecMesonCut->passTrackEtaEast(lTrackA)) // trackB
          //  {
          //    Float_t  w = mVecMesonCorr->getWeight(lTrackA);
          //    TVector2 q2VectorA = mVecMesonCorr->calq2Vector(lTrackA);
          //    TVector2 q2CorrA   = mVecMesonCorr->getReCenterPar_East(mCent,mRunIndex,mVzSign);
          //    Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);
          //  }
          //  if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
          //  {
          //    Float_t  w = mVecMesonCorr->getWeight(lTrackB);
          //    TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
          //    TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(mCent,mRunIndex,mVzSign);
          //    Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
          //  }
          //  Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
          //  Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,mRunIndex,mCent,mVzSign);
          //  Float_t PhiMinusPsi = lTrack.Phi() - Psi2_east;

          //  if(PhiMinusPsi >= -(3.0/2.0)*TMath::Pi() && PhiMinusPsi < -(1.0/2.0)*TMath::Pi()) PhiMinusPsi += TMath::Pi();
          //  if(PhiMinusPsi <=  (3.0/2.0)*TMath::Pi() && PhiMinusPsi >  (1.0/2.0)*TMath::Pi()) PhiMinusPsi -= TMath::Pi();
          //  PhiMinusPsi = TMath::Abs(PhiMinusPsi);

          //  mVecMesonHistoFlow->FillSys(pt_lTrack,mCent,PhiMinusPsi,i_dca,i_sig,i_nhit,Res2,InvMass_lTrack,reweight,mX_flag,2);
          //}          
          //////////////// ELLIPTIC FLOW ///////////////////////////////////////////////////                
        //}
      }
    }
  }
}


void StVecMesonTree::FillKStar(StMesonTrack* Meson_track, Int_t mCent, Int_t vz, Int_t Psi2)
{

  // Initialise Track 
  Float_t m2A = -999.9;
  Float_t m2B = -999.9;
  Float_t nsA = -999.9;
  Float_t nsB = -999.9;
  Float_t dcaA = -999.9;
  Float_t dcaB = -999.9;
  Float_t nhA = -1;
  Float_t nhB = -1;
  TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
  TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
  Int_t flagA = -1;
  Int_t flagB = -1;

  TVector2 mQv2West = mQ2West[mCent][vz][Psi2][0];
  TVector2 mQv2East = mQ2East[mCent][vz][Psi2][0];

  int mVzSign = mVzEPBin[mCent][vz][Psi2][0];
  int mRunIndex = mRunIdx[mCent][vz][Psi2][0];
  int reweight = mReweight[mCent][vz][Psi2][0];

  // get Track Information
  if(mVecMesonCorr->passTrackEtaNumCut(mNumTrackEast[mCent][vz][Psi2][0],mNumTrackWest[mCent][vz][Psi2][0]))
  {
    m2A = Meson_track->getMass2A();
    m2B = Meson_track->getMass2B();
    nsA = Meson_track->getNSigA();
    nsB = Meson_track->getNSigB();
    dcaA = Meson_track->getDcaA();
    dcaB = Meson_track->getDcaB();
    nhA = Meson_track->getNHitA();
    nhB = Meson_track->getNHitB();
    lTrackA = Meson_track->getTrackA();
    lTrackB = Meson_track->getTrackB();
    flagA = Meson_track->getFlagA();
    flagB = Meson_track->getFlagB();
    Float_t pA = lTrackA.P();
    Float_t pB = lTrackB.P();
    TLorentzVector lTrack = lTrackA + lTrackB; // KStar-meson
    Float_t pt_lTrack = lTrack.Perp();
    Float_t rapidity_lTrack = lTrack.Rapidity();

    if(TMath::Abs(rapidity_lTrack) > vmsa::mEtaMax) return;

    Float_t InvMass_lTrack = lTrack.M();
    TVector3 vBetaKStar = -1.0*lTrack.BoostVector(); // get phi beta
    TLorentzVector lKRest = lTrackA;
    lKRest.Boost(vBetaKStar); // boost kaon back to phi rest frame
    TVector3 vKRest = lKRest.Vect().Unit(); // kaon momentum direction in phi rest frame
    
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++) // systematic loop for dca
    {
      if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,2)) ) continue;
      mVecMesonHistoManger->FillDcaSys(dcaA,dcaB,i_dca); // fill QA for dcaA and dcaB

      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++) // systematic loop for nSigma
      {
        if(i_dca != 0 && i_sig != 0) continue;
        if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,2)) ) continue;
        mVecMesonHistoManger->FillSigSys(nsA,nsB,i_sig); // fill QA for nsA and nsB 
    
        for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++) // systematic loop for nSigma
        {
          if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
          if( !(mVecMesonCut->passTrackNHitSys(nhA,nhB,i_nhit,2)) ) continue;
          mVecMesonHistoManger->FillNHitSys(nhA,nhB,i_nhit); // fill QA for nhA and nhB 

          //////////////// GLOBAL SPIN ALIGNMENT ///////////////////////////////////////////////////
          if(mVecMesonCut->passEtaEast(lTrackA)) // K+/- neg eta(east)
          { // Below is West Only
            TVector2 Q2Vector = mQv2West;
            // subtract auto-correlation from pos eta(west) event plane
            if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackB);
              TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
              TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
            }
            Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
            Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,mRunIndex,mCent,mVzSign);

            TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
            TVector3 nQ = nQ_West.Unit();
            Double_t CosThetaStar = vKRest.Dot(nQ);
            mVecMesonHistoManger->FillSys(pt_lTrack,mCent,CosThetaStar,i_dca,i_sig,i_nhit,Res2,InvMass_lTrack,reweight,mX_flag,2);

            //TVector3 nQ_West_EP(TMath::Cos(Psi2_west),TMath::Sin(Psi2_west),0.0); // tangent vector of 2nd Event Plane
            //TVector3 nQ_EP = nQ_West_EP.Unit();          
            //Double_t CosThetaStar_EP = vKRest.Dot(nQ_EP);
          }
    
          if(mVecMesonCut->passEtaWest(lTrackA)) // K+/- pos eta (west)
          { // Below is East Only
            TVector2 Q2Vector = mQv2East;
            // subtract auto-correlation from pos eta(west) event plane
            if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackB);
              TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
              TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
            }
            Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
            Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,mRunIndex,mCent,mVzSign);
            
            TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane
            TVector3 nQ = nQ_East.Unit();
            Double_t CosThetaStar = vKRest.Dot(nQ);
            mVecMesonHistoManger->FillSys(pt_lTrack,mCent,CosThetaStar,i_dca,i_sig,i_nhit,Res2,InvMass_lTrack,reweight,mX_flag,2);
            
            //TVector3 nQ_East_EP(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0); // tangent vector of 2nd Event Plane
            //TVector3 nQ_EP = nQ_East_EP.Unit();
            //Double_t CosThetaStar_EP = vKRest.Dot(nQ_EP);
            //mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,2);
          }
          //////////////// GLOBAL SPIN ALIGNMENT ///////////////////////////////////////////////////
          
          //////////////// ELLIPTIC FLOW ///////////////////////////////////////////////////
          if(mVecMesonCut->passEtaEast(lTrack)) // phi neg eta(east)
          { // Below is West Only
            TVector2 Q2Vector = mQv2West;
            // subtract auto-correlation from pos eta(west) event plane
            if(flagA == 0 && mVecMesonCut->passTrackEP(lTrackA,dcaA) && mVecMesonCut->passTrackEtaWest(lTrackA)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackA);
              TVector2 q2VectorA = mVecMesonCorr->calq2Vector(lTrackA);
              TVector2 q2CorrA   = mVecMesonCorr->getReCenterPar_West(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);
            }
            if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackB);
              TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
              TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
            }
            Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
            Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,mRunIndex,mCent,mVzSign);
            Float_t PhiMinusPsi = lTrack.Phi() - Psi2_west;
            
            if(PhiMinusPsi >= -(3.0/2.0)*TMath::Pi() && PhiMinusPsi < -(1.0/2.0)*TMath::Pi()) PhiMinusPsi += TMath::Pi();
            if(PhiMinusPsi <=  (3.0/2.0)*TMath::Pi() && PhiMinusPsi >  (1.0/2.0)*TMath::Pi()) PhiMinusPsi -= TMath::Pi();
            PhiMinusPsi = TMath::Abs(PhiMinusPsi);
           
            mVecMesonHistoFlow->FillSys(pt_lTrack,mCent,PhiMinusPsi,i_dca,i_sig,i_nhit,Res2,InvMass_lTrack,reweight,mX_flag,2);
          }
  
          if(mVecMesonCut->passEtaWest(lTrack)) // phi pos eta (west)
          { // Below is East Only
            TVector2 Q2Vector = mQv2East;
            // subtract auto-correlation from pos eta(west) event plane
            if(flagA == 0 && mVecMesonCut->passTrackEP(lTrackA,dcaA) && mVecMesonCut->passTrackEtaEast(lTrackA)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackA);
              TVector2 q2VectorA = mVecMesonCorr->calq2Vector(lTrackA);
              TVector2 q2CorrA   = mVecMesonCorr->getReCenterPar_East(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);
            }
            if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
            {
              Float_t  w = mVecMesonCorr->getWeight(lTrackB);
              TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
              TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(mCent,mRunIndex,mVzSign);
              Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
            }
            Float_t Res2 = mVecMesonCorr->getResolution2_EP(mCent);
            Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,mRunIndex,mCent,mVzSign);
            Float_t PhiMinusPsi = lTrack.Phi() - Psi2_east;

            if(PhiMinusPsi >= -(3.0/2.0)*TMath::Pi() && PhiMinusPsi < -(1.0/2.0)*TMath::Pi()) PhiMinusPsi += TMath::Pi();
            if(PhiMinusPsi <=  (3.0/2.0)*TMath::Pi() && PhiMinusPsi >  (1.0/2.0)*TMath::Pi()) PhiMinusPsi -= TMath::Pi();
            PhiMinusPsi = TMath::Abs(PhiMinusPsi);

            mVecMesonHistoFlow->FillSys(pt_lTrack,mCent,PhiMinusPsi,i_dca,i_sig,i_nhit,Res2,InvMass_lTrack,reweight,mX_flag,2);
          }          
          //////////////// ELLIPTIC FLOW ///////////////////////////////////////////////////                
        }
      }
    }
  }
}

	    
//------------------------------------------------------------------------------------------------------------------

/*void StVecMesonTree::MixEvent_Rho(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2)
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
  //mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  //mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  //mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  //mEventId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->eventId()));
  //mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  //mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  //mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  //mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  //mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  //mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  //mQ2East[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2East));
  //mQ2West[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2West));
  //mQ2Full[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2Full));
  //mNumTrackEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaEast));
  //mNumTrackWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaWest));
  //mNumTrackFull[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFull));
  //mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullEast));
  //mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullWest));

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    //if(mVecMesonCut->passTrackMeson(track,event))
    {
      Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track,pico);
      Float_t scale_nSigma_factor = vmsa::mSigScaleMap[mEnergy];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().Mag();
      Float_t Mass2_low = -0.12;
      Float_t Mass2_up = 0.16;

      Int_t charge = 0; // pi+
      if(Polarity < 0) charge = 1; // pi-

      if(mVecMesonCut->passSigPionCut(track,scale_nSigma_factor))
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
	  mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaPion())*scale_nSigma_factor)); // nSigma
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
*/
//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // start to select phi candidate in a event
     // cout << "Entered Loop" << endl;    
 
      MEKey key_kplus   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0); //K+
      MEKey key_kminus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1); //K-

      TLorentzVector ltrackA, ltrackB;
      for(Int_t n_kplus = 0; n_kplus < mHelix[key_kplus].size(); n_kplus++) // first track loop over K+ candidates
      {
       // cout << "Loop over K+" << endl;
	StThreeVectorF p_vecA = mHelix[key_kplus][n_kplus].cat(mHelix[key_kplus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	p_vecA *= mMomentum[key_kplus][n_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_kminus].size(); n_kminus++) // second track loop over K- candidates
	{
         // cout << "Loop over K-" << endl;
	  StThreeVectorF p_vecB = mHelix[key_kminus][n_kminus].cat(mHelix[key_kminus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_kminus][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);
	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

          
          //cout << "About to fill phi Meson" << endl;
	  if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	  {
            StMesonTrack *mesonTrack = new StMesonTrack(); 
	    //mesonTrack->setMass2A(mMass2[key_kplus][n_kplus]); // K+
	    //mesonTrack->setMass2B(mMass2[key_kminus][n_kminus]); // K-
	    mesonTrack->setNSigA(mNSigma[key_kplus][n_kplus]); // K+
	    mesonTrack->setNSigB(mNSigma[key_kminus][n_kminus]); // K-
	    mesonTrack->setDcaA(mDca[key_kplus][n_kplus]); // K+
	    mesonTrack->setDcaB(mDca[key_kminus][n_kminus]); // K-
            //mesonTrack->setNHitA(mNHitsFit[key_kplus][n_kplus]); // K+
            //mesonTrack->setNHitB(mNHitsFit[key_kminus][n_kminus]); // K-
	    mesonTrack->setTrackA(ltrackA); // K+
	    mesonTrack->setTrackB(ltrackB); // K-
	    mesonTrack->setFlagA(Bin_Event); // K+
	    mesonTrack->setFlagB(Bin_Event); // K-
	    h_Mass2->Fill(pt,InvMassAB);
            FillPhi(mesonTrack,cent9,Bin_vz,Bin_Psi2);
            delete mesonTrack;
	  }
	}
      }
    }
  }

  if(Flag_ME == 1) // mixed event
  {
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_kplus   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0); //K+
      MEKey key_A_kminus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1); //K-
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
        MEKey key_B_kplus   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0); //K+
        MEKey key_B_kminus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1); //K-


	TLorentzVector ltrackA, ltrackB;

	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	for(Int_t n_kplus = 0; n_kplus < mHelix[key_A_kplus].size(); n_kplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_kplus][n_kplus].cat(mHelix[key_A_kplus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_kplus][n_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+

	  for(Int_t n_kminus = 0; n_kminus < mHelix[key_B_kminus].size(); n_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_kminus][n_kminus].cat(mHelix[key_B_kminus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_kminus][n_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
              StMesonTrack *mesonTrack = new StMesonTrack(); 
	      //mesonTrack->setMass2A(mMass2[key_A_kplus][n_kplus]); // K+
	      //mesonTrack->setMass2B(mMass2[key_B_kminus][n_kminus]); // K-
	      mesonTrack->setNSigA(mNSigma[key_A_kplus][n_kplus]); // K+
	      mesonTrack->setNSigB(mNSigma[key_B_kminus][n_kminus]); // K-
	      mesonTrack->setDcaA(mDca[key_A_kplus][n_kplus]); // K+
	      mesonTrack->setDcaB(mDca[key_B_kminus][n_kminus]); // K-
              //mesonTrack->setNHitA(mNHitsFit[key_A_kplus][n_kplus]); // K+
              //mesonTrack->setNHitB(mNHitsFit[key_B_kminus][n_kminus]); // K-
	      mesonTrack->setTrackA(ltrackA); // K+
	      mesonTrack->setTrackB(ltrackB); // K-
	      mesonTrack->setFlagA(Bin_Event_A); // K+
	      mesonTrack->setFlagB(Bin_Event_B); // K-
	      h_Mass2->Fill(pt,InvMassAB);
              FillPhi(mesonTrack,cent9,Bin_vz,Bin_Psi2);
              delete mesonTrack;
	    }
	  }
	}

        // mix pi- candidates from A event with K+ candidates from B event
	for(Int_t n_kminus = 0; n_kminus < mHelix[key_A_kminus].size(); n_kminus++) // first track loop over pi- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix[key_A_kminus][n_kminus].cat(mHelix[key_A_kminus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_kminus][n_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-

	  for(Int_t n_kplus = 0; n_kplus < mHelix[key_B_kplus].size(); n_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix[key_B_kplus][n_kplus].cat(mHelix[key_B_kplus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_kplus][n_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
              StMesonTrack *mesonTrack = new StMesonTrack(); 
	      //mesonTrack->setMass2A(mMass2[key_B_kplus][n_kplus]); // K+
	      //mesonTrack->setMass2B(mMass2[key_A_kminus][n_kminus]); // K-
	      mesonTrack->setNSigA(mNSigma[key_B_kplus][n_kplus]); // K+
	      mesonTrack->setNSigB(mNSigma[key_A_kminus][n_kminus]); // K-
	      mesonTrack->setDcaA(mDca[key_B_kplus][n_kplus]); // K+
	      mesonTrack->setDcaB(mDca[key_A_kminus][n_kminus]); // K-
              //mesonTrack->setNHitA(mNHitsFit[key_B_kplus][n_kplus]); // K+
              //mesonTrack->setNHitB(mNHitsFit[key_A_kminus][n_kminus]); // K-
	      mesonTrack->setTrackA(ltrackB); // K+
	      mesonTrack->setTrackB(ltrackA); // K-
	      mesonTrack->setFlagA(Bin_Event_B); // K+
	      mesonTrack->setFlagB(Bin_Event_A); // K-
	      h_Mass2->Fill(pt,InvMassAB);
              FillPhi(mesonTrack,cent9,Bin_vz,Bin_Psi2);
              delete mesonTrack;
	    }
	  }
	}
      }
    }
  }
}

void StVecMesonTree::doKStar(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
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
            StMesonTrack *mesonTrack = new StMesonTrack(); 
	    mesonTrack->setMass2A(mMass2[key_kplus][n_kplus]); // K+
	    mesonTrack->setMass2B(mMass2[key_piminus][n_piminus]); // pi-
	    mesonTrack->setNSigA(mNSigma[key_kplus][n_kplus]); // K+
	    mesonTrack->setNSigB(mNSigma[key_piminus][n_piminus]); // pi-
	    mesonTrack->setDcaA(mDca[key_kplus][n_kplus]); // K+
	    mesonTrack->setDcaB(mDca[key_piminus][n_piminus]); // pi-
            mesonTrack->setNHitA(mNHitsFit[key_kplus][n_kplus]); // K+
            mesonTrack->setNHitB(mNHitsFit[key_piminus][n_piminus]); // pi-
	    mesonTrack->setTrackA(ltrackA); // K+
	    mesonTrack->setTrackB(ltrackB); // pi-
	    mesonTrack->setFlagA(Bin_Event); // K+
	    mesonTrack->setFlagB(Bin_Event); // pi-
	    h_Mass2->Fill(pt,InvMassAB);
            FillKStar(mesonTrack,cent9,Bin_vz,Bin_Psi2);
            delete mesonTrack;
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
            StMesonTrack *mesonTrack = new StMesonTrack(); 
	    mesonTrack->setMass2B(mMass2[key_piplus][n_piplus]); // pi+
	    mesonTrack->setMass2A(mMass2[key_kminus][n_kminus]); // K-
	    mesonTrack->setNSigB(mNSigma[key_piplus][n_piplus]); // pi+
	    mesonTrack->setNSigA(mNSigma[key_kminus][n_kminus]); // K-
	    mesonTrack->setDcaB(mDca[key_piplus][n_piplus]); // pi+
	    mesonTrack->setDcaA(mDca[key_kminus][n_kminus]); // K-
            mesonTrack->setNHitB(mNHitsFit[key_piplus][n_piplus]); // pi+
            mesonTrack->setNHitA(mNHitsFit[key_kminus][n_kminus]); // K-
	    mesonTrack->setTrackB(ltrackA); // pi+
	    mesonTrack->setTrackA(ltrackB); // K-
	    mesonTrack->setFlagB(Bin_Event); // pi+
	    mesonTrack->setFlagA(Bin_Event); // K-
	    h_Mass2->Fill(pt,InvMassAB);
            FillKStar(mesonTrack,cent9,Bin_vz,Bin_Psi2);
            delete mesonTrack;
	  }
	}
      }
    }
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
              StMesonTrack *mesonTrack = new StMesonTrack(); 
	      mesonTrack->setMass2A(mMass2[key_A_kplus][n_kplus]); // K+
	      mesonTrack->setMass2B(mMass2[key_B_piminus][n_piminus]); // pi-
	      mesonTrack->setNSigA(mNSigma[key_A_kplus][n_kplus]); // K+
	      mesonTrack->setNSigB(mNSigma[key_B_piminus][n_piminus]); // pi-
	      mesonTrack->setDcaA(mDca[key_A_kplus][n_kplus]); // K+
	      mesonTrack->setDcaB(mDca[key_B_piminus][n_piminus]); // pi-
              mesonTrack->setNHitA(mNHitsFit[key_A_kplus][n_kplus]); // K+
              mesonTrack->setNHitB(mNHitsFit[key_B_piminus][n_piminus]); // pi-
	      mesonTrack->setTrackA(ltrackA); // K+
	      mesonTrack->setTrackB(ltrackB); // pi-
	      mesonTrack->setFlagA(Bin_Event_A); // K+
	      mesonTrack->setFlagB(Bin_Event_B); // pi-
	      h_Mass2->Fill(pt,InvMassAB);
              FillKStar(mesonTrack,cent9,Bin_vz,Bin_Psi2);
              delete mesonTrack;
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
              StMesonTrack *mesonTrack = new StMesonTrack(); 
	      mesonTrack->setMass2B(mMass2[key_A_piplus][n_piplus]); // pi+
	      mesonTrack->setMass2A(mMass2[key_B_kminus][n_kminus]); // K-
	      mesonTrack->setNSigB(mNSigma[key_A_piplus][n_piplus]); // pi+
	      mesonTrack->setNSigA(mNSigma[key_B_kminus][n_kminus]); // K-
	      mesonTrack->setDcaB(mDca[key_A_piplus][n_piplus]); // pi+
	      mesonTrack->setDcaA(mDca[key_B_kminus][n_kminus]); // K-
              mesonTrack->setNHitB(mNHitsFit[key_A_piplus][n_piplus]); // pi+
              mesonTrack->setNHitA(mNHitsFit[key_B_kminus][n_kminus]); // K-
	      mesonTrack->setTrackB(ltrackA); // pi+
	      mesonTrack->setTrackA(ltrackB); // K-
	      mesonTrack->setFlagB(Bin_Event_A); // pi+
	      mesonTrack->setFlagA(Bin_Event_B); // K-
	      h_Mass2->Fill(pt,InvMassAB);
              FillKStar(mesonTrack,cent9,Bin_vz,Bin_Psi2);
              delete mesonTrack;
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
              StMesonTrack *mesonTrack = new StMesonTrack(); 
	      mesonTrack->setMass2B(mMass2[key_B_piplus][n_piplus]); // pi+
	      mesonTrack->setMass2A(mMass2[key_A_kminus][n_kminus]); // K-
	      mesonTrack->setNSigB(mNSigma[key_B_piplus][n_piplus]); // pi+
	      mesonTrack->setNSigA(mNSigma[key_A_kminus][n_kminus]); // K-
	      mesonTrack->setDcaB(mDca[key_B_piplus][n_piplus]); // pi+
	      mesonTrack->setDcaA(mDca[key_A_kminus][n_kminus]); // K-
              mesonTrack->setNHitB(mNHitsFit[key_B_piplus][n_piplus]); // pi+
              mesonTrack->setNHitA(mNHitsFit[key_A_kminus][n_kminus]); // K-
	      mesonTrack->setTrackB(ltrackB); // pi+
	      mesonTrack->setTrackA(ltrackA); // K-
	      mesonTrack->setFlagB(Bin_Event_B); // pi+
	      mesonTrack->setFlagA(Bin_Event_A); // K-
	      h_Mass2->Fill(pt,InvMassAB);
              FillKStar(mesonTrack,cent9,Bin_vz,Bin_Psi2);
              delete mesonTrack;
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
              StMesonTrack *mesonTrack = new StMesonTrack(); 
	      mesonTrack->setMass2A(mMass2[key_B_kplus][n_kplus]); // K+
	      mesonTrack->setMass2B(mMass2[key_A_piminus][n_piminus]); // pi-
	      mesonTrack->setNSigA(mNSigma[key_B_kplus][n_kplus]); // K+
	      mesonTrack->setNSigB(mNSigma[key_A_piminus][n_piminus]); // pi-
	      mesonTrack->setDcaA(mDca[key_B_kplus][n_kplus]); // K+
	      mesonTrack->setDcaB(mDca[key_A_piminus][n_piminus]); // pi-
              mesonTrack->setNHitA(mNHitsFit[key_B_kplus][n_kplus]); // K+
              mesonTrack->setNHitB(mNHitsFit[key_A_piminus][n_piminus]); // pi-
	      mesonTrack->setTrackA(ltrackB); // K+
	      mesonTrack->setTrackB(ltrackA); // pi-
	      mesonTrack->setFlagA(Bin_Event_B); // K+
	      mesonTrack->setFlagB(Bin_Event_A); // pi-
	      h_Mass2->Fill(pt,InvMassAB);
              FillKStar(mesonTrack,cent9,Bin_vz,Bin_Psi2);
              delete mesonTrack;
	    }
	  }
	}
      }
    }
  }
}


void StVecMesonTree::MixEvent_Phi(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2, Int_t runIndex, Float_t reweight, Int_t vzEPbin)
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
  const Int_t nTracks = pico->numberOfTracks();

  // store Event Information
  StThreeVectorF primVer; 
  primVer.setX(event->primaryVertex().x());
  primVer.setY(event->primaryVertex().y());
  primVer.setZ(event->primaryVertex().z());

  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(primVer));
  mRunIdx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(runIndex));
  mVzEPBin[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(vzEPbin));
  mReweight[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(reweight));
  //mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  //mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  //mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  //mEventId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->eventId()));
  //mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  //mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  //mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  //mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  //mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  //mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  mQ2East[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2East));
  mQ2West[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2West));
  //mQ2Full[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2Full));
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaEast));
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaWest));
  //mNumTrackFull[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFull));
  //mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullEast));
  //mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullWest));

  //cout << "filled event information" << endl;
  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mVecMesonCut->passTrackMeson(track,event,2))
    {
      //Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track,pico);
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().Mag();
      //Float_t Mass2K_low = 0.16;
      //Float_t Mass2K_up = 0.36;
      //Float_t Mass2Pi_low = -0.20;
      //Float_t Mass2Pi_up = 0.15;
      
      Int_t charge = 0; // +
      if(Polarity < 0) charge = 1; // -
     
      //cout << "grabbed some track info" << endl;
      Float_t nSigmaKaon = track->nSigmaKaon();
      //cout << "grabbed nSigmaKaon = " << nSigmaKaon << endl;
      {
	if(
            (fabs(nSigmaKaon) < 3.0)
            //(Mass2 > Mass2K_low && Mass2 < Mass2K_up) && fabs(nSigmaKaon) < 2.0
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
          //cout << "about to create key" << endl;
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
          //cout << "create key" << endl;
	  //mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
          //cout << "pushback mass" << endl;
	  mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
          //cout << "pushback mass" << endl;
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
          //cout << "pushback mass" << endl;
	  mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaKaon()))); // nSigma
          //cout << "pushback mass" << endl;
          //key.print();
    	  mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
          //cout << "pushback mass" << endl;
	  mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
          //cout << "pushback mass" << endl;
          //mCharge[key].push_back(static_cast<Int_t>(charge));
          h_KdEdxRig->Fill(momentum*Polarity,track->dEdx());
          //cout << "fill" << endl;
          //h_KM2Rig->Fill(momentum*Polarity,Mass2);
          //cout << "fill" << endl;
          //h_KInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
          //cout << "something too large in kaon container" << endl;
          //cout << "push back kaon info" << endl;
	}
        else continue;
      }
      //Float_t nSigmaPion = track->nSigmaPion();
      //{
      //  charge += 2; // changes id to pi+ and pi-
      //  if(
      //      (Mass2 > Mass2Pi_low && Mass2 < Mass2Pi_up) || (Mass2 < -10.0 && fabs(nSigmaPion) < 2.0) 
      //      //(Mass2 > Mass2Pi_low && Mass2 < Mass2Pi_up) && fabs(nSigmaPion) < 2.0
      //    )
      //  {
      //    StThreeVectorD primMom; 
      //    primMom.setX(track->pMom().x());
      //    primMom.setY(track->pMom().y());
      //    primMom.setZ(track->pMom().z());
      //    StThreeVectorD primVer; 
      //    primVer.setX(event->primaryVertex().x());
      //    primVer.setY(event->primaryVertex().y());
      //    primVer.setZ(event->primaryVertex().z());
      //    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
      //    mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
      //    mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
      //    mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
      //    mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaPion()))); // nSigma
      //    mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
      //    mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
      //    //mCharge[key].push_back(static_cast<Int_t>(charge));
      //    
      //    h_PidEdxRig->Fill(momentum*Polarity,track->dEdx());
      //    h_PiM2Rig->Fill(momentum*Polarity,Mass2);
      //    h_PiInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
      //    cout << "push back pion info" << endl;
      //    //cout << "something too large in pion contianer" << endl;
      //  }
      //}

    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_phi(cent9,Bin_vz,Bin_Psi2,0);
  }

  if(Flag_ME == 1) // mix event
  {
    //cout << "about to mix event" << endl;
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth)
    {
      doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      //cout << "mix the event" << endl;
      clear_phi(cent9,Bin_vz,Bin_Psi2,1);
    }
  }

}


void StVecMesonTree::MixEvent_KStar(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2, Int_t runIndex, Float_t reweight, Int_t vzEPbin)
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

  // store Event Information
  StThreeVectorF primVer; 
  primVer.setX(event->primaryVertex().x());
  primVer.setY(event->primaryVertex().y());
  primVer.setZ(event->primaryVertex().z());

  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(primVer));
  mRunIdx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(runIndex));
  mVzEPBin[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(vzEPbin));
  mReweight[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(reweight));
  //mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  //mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  //mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  //mEventId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->eventId()));
  //mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  //mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  //mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  //mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  //mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  //mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  mQ2East[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2East));
  mQ2West[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2West));
  //mQ2Full[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2Full));
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaEast));
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaWest));
  //mNumTrackFull[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFull));
  //mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullEast));
  //mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullWest));

  cout << "filled event information" << endl;
  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mVecMesonCut->passTrackMeson(track,event,2))
    {
      Float_t Mass2 = mVecMesonCut->getPrimaryMass2(track,pico);
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().Mag();
      Float_t Mass2K_low = 0.16;
      Float_t Mass2K_up = 0.36;
      Float_t Mass2Pi_low = -0.20;
      Float_t Mass2Pi_up = 0.15;
      
      Int_t charge = 0; // +
      if(Polarity < 0) charge = 1; // -
     
      cout << "grabbed some track info" << endl;
      Float_t nSigmaKaon = track->nSigmaKaon();
      cout << "grabbed nSigmaKaon = " << nSigmaKaon << endl;
      {
	if(
            (Mass2 > Mass2K_low && Mass2 < Mass2K_up) || (Mass2 < -10.0 && fabs(nSigmaKaon) < 2.0)
            //(Mass2 > Mass2K_low && Mass2 < Mass2K_up) && fabs(nSigmaKaon) < 2.0
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
          cout << "about to create key" << endl;
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
          cout << "create key" << endl;
	  mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
          cout << "pushback mass" << endl;
	  mDca[key].push_back(static_cast<Float_t>(track->gDCA(primVer.x(),primVer.y(),primVer.z())*Polarity)); // dca*charge 
          cout << "pushback mass" << endl;
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
          cout << "pushback mass" << endl;
	  mNSigma[key].push_back(static_cast<Float_t>((track->nSigmaKaon()))); // nSigma
          cout << "pushback mass" << endl;
          //key.print();
    	  mHelix[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(primMom,primVer,event->bField()*MAGFIELDFACTOR,Polarity)));// get helix from the pMom 
          cout << "pushback mass" << endl;
	  mMomentum[key].push_back(static_cast<Float_t>(momentum));// get helix from the pMom 
          cout << "pushback mass" << endl;
          //mCharge[key].push_back(static_cast<Int_t>(charge));
          h_KdEdxRig->Fill(momentum*Polarity,track->dEdx());
          cout << "fill" << endl;
          h_KM2Rig->Fill(momentum*Polarity,Mass2);
          cout << "fill" << endl;
          h_KInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
          //cout << "something too large in kaon container" << endl;
          cout << "push back kaon info" << endl;
	}
      else continue;
      }
      Float_t nSigmaPion = track->nSigmaPion();
      {
        charge += 2; // changes id to pi+ and pi-
	if(
            (Mass2 > Mass2Pi_low && Mass2 < Mass2Pi_up) || (Mass2 < -10.0 && fabs(nSigmaPion) < 2.0) 
            //(Mass2 > Mass2Pi_low && Mass2 < Mass2Pi_up) && fabs(nSigmaPion) < 2.0
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
          //mCharge[key].push_back(static_cast<Int_t>(charge));
          
          h_PidEdxRig->Fill(momentum*Polarity,track->dEdx());
          h_PiM2Rig->Fill(momentum*Polarity,Mass2);
          h_PiInvBetaRig->Fill(momentum*Polarity,1.0/mVecMesonCut->getBeta(track,pico));
          cout << "push back pion info" << endl;
          //cout << "something too large in pion contianer" << endl;
	}
      }

    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doKStar(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_kstar(cent9,Bin_vz,Bin_Psi2,0);
  }

  if(Flag_ME == 1) // mix event
  {
    cout << "about to mix event" << endl;
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth)
    {
      doKStar(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      cout << "mix the event" << endl;
      clear_kstar(cent9,Bin_vz,Bin_Psi2,1);
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
 
  //mRunIndex = -1;
  //mReweight = 1.0;
  //mVzSign = -1;
  //mCent = -1;
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

//void StVecMesonTree::passExtraEvent(Int_t RunIndex, Int_t Centrality, Int_t VzSign, Double_t Reweight)
//{
//  mRunIndex = RunIndex;
//  mCent = Centrality;
//  mVzSign = VzSign;
//  mReweight = Reweight;
//}
//------------------------------------------------------------------------------------------------------------------
