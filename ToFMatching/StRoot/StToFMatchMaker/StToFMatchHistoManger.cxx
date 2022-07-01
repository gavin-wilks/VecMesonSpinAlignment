#include "StRoot/StToFMatchMaker/StToFMatchHistoManger.h"
#include "StRoot/StToFMatchMaker/StToFMatchCons.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include <string>

ClassImp(StToFMatchHistoManger)

//-------------------------------------------------------------------------------------------

StToFMatchHistoManger::StToFMatchHistoManger()
{
}

//-------------------------------------------------------------------------------------------

StToFMatchHistoManger::~StToFMatchHistoManger()
{
  /* */
}

//-------------------------------------------------------------------------------------------
void StToFMatchHistoManger::InitQA()
{
  h_mVz   = new TH1F("h_mVz","h_mVz",201,-100.5,100.5);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",1000,-0.5,999.5);

  h_mDEdx = new TH2F("h_mDEdx","h_mDEdx",1000,-10.0,10.0,1000,0,40);
  h_mMass2 = new TH2F("h_mMass2","h_mMass2",1000,-10.0,10.0,1000,-0.3,1.7);

  h_mDEdx_Pion = new TH2F("h_mDEdx_Pion","h_mDEdx_Pion",1000,-10.0,10.0,1000,0,40);
  h_mMass2_Pion = new TH2F("h_mMass2_Pion","h_mMass2_Pion",1000,-10.0,10.0,1000,-0.3,1.7);

  h_mDEdx_Kaon = new TH2F("h_mDEdx_Kaon","h_mDEdx_Kaon",1000,-10.0,10.0,1000,0,40);
  h_mMass2_Kaon = new TH2F("h_mMass2_Kaon","h_mMass2_Kaon",1000,-10.0,10.0,1000,-0.3,1.7);

  h_mDEdx_Proton = new TH2F("h_mDEdx_Proton","h_mDEdx_Proton",1000,-10.0,10.0,1000,0,40);
  h_mMass2_Proton = new TH2F("h_mMass2_Proton","h_mMass2_Proton",1000,-10.0,10.0,1000,-0.3,1.7);
}

void StToFMatchHistoManger::FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx->Fill(p,dEdx);
  h_mMass2->Fill(p,Mass2);
}

void StToFMatchHistoManger::FillQA_Pion(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx_Pion->Fill(p,dEdx);
  h_mMass2_Pion->Fill(p,Mass2);
}

void StToFMatchHistoManger::FillQA_Kaon(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx_Kaon->Fill(p,dEdx);
  h_mMass2_Kaon->Fill(p,Mass2);
}

void StToFMatchHistoManger::FillQA_Proton(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx_Proton->Fill(p,dEdx);
  h_mMass2_Proton->Fill(p,Mass2);
}

void StToFMatchHistoManger::FillQA_Event(Float_t vz, Float_t refMult)
{
  h_mVz->Fill(vz);
  h_mRefMult->Fill(refMult);
}

void StToFMatchHistoManger::WriteQA()
{
  h_mVz->Write();
  h_mRefMult->Write();

  h_mDEdx->Write();
  h_mMass2->Write();

  h_mDEdx_Pion->Write();
  h_mMass2_Pion->Write();

  h_mDEdx_Kaon->Write();
  h_mMass2_Kaon->Write();

  // h_mDEdx_Proton->Write();
  // h_mMass2_Proton->Write();
}
//-------------------------------------------------------------------------------------------
void StToFMatchHistoManger::InitHist()
{
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_TPC[i_charge][i_pid][i_cent] = new TH3D(HistName_Trakcs_TPC.c_str(),HistName_Trakcs_TPC.c_str(),tof::BinPt,tof::ptMin,tof::ptMax,tof::BinEta,tof::etaMin,tof::etaMax,tof::BinPhi,tof::phiMin,tof::phiMax);
	h_mTracks_TPC[i_charge][i_pid][i_cent]->Sumw2();

	string HistName_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_ToF[i_charge][i_pid][i_cent] = new TH3D(HistName_Trakcs_ToF.c_str(),HistName_Trakcs_ToF.c_str(),tof::BinPt,tof::ptMin,tof::ptMax,tof::BinEta,tof::etaMin,tof::etaMax,tof::BinPhi,tof::phiMin,tof::phiMax);
	h_mTracks_ToF[i_charge][i_pid][i_cent]->Sumw2();

	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_TPC[HistName_TPC] = new TH1D(HistName_TPC.c_str(),HistName_TPC.c_str(),tof::BinPt,tof::ptMin,tof::ptMax);
	    h_mCounts_TPC[HistName_TPC]->Sumw2();
	    string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_ToF[HistName_ToF] = new TH1D(HistName_ToF.c_str(),HistName_ToF.c_str(),tof::BinPt,tof::ptMin,tof::ptMax);
	    h_mCounts_ToF[HistName_ToF]->Sumw2();
	  }
	}
      }
    }
  }
  h_FrameEta_ToF = new TH1D("h_FrameEta_ToF","h_FrameEta_ToF",tof::BinEta,tof::etaMin,tof::etaMax);
  h_FramePhi_ToF = new TH1D("h_FramePhi_ToF","h_FramePhi_ToF",tof::BinPhi,tof::phiMin,tof::phiMax);
}

int StToFMatchHistoManger::findEta(float eta)
{
  int EtaBin = h_FrameEta_ToF->FindBin(eta)-1;
  return EtaBin;
}

int StToFMatchHistoManger::findPhi(float phi)
{
  int PhiBin = h_FramePhi_ToF->FindBin(phi)-1;
  return PhiBin;
}

void StToFMatchHistoManger::Fill_TPC(int charge, int pid, int cent, float pt, float eta, float phi)
{
  h_mTracks_TPC[charge][pid][cent]->Fill(pt,eta,phi);
  h_mTracks_TPC[charge][pid][9]->Fill(pt,eta,phi);  //miniBias
  int i_eta = findEta(eta);
  int i_phi = findPhi(phi);
  string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
  h_mCounts_TPC[HistName_TPC]->Fill(pt);

  string HistName_TPC_miniBias = Form("h_mCounts_TPC_%s%s_Cent_9_Eta_%d_Phi_%d",tof::mPID_ToF[pid].c_str(),tof::mCharge[charge].c_str(),i_eta,i_phi);
  h_mCounts_TPC[HistName_TPC_miniBias]->Fill(pt);
}

void StToFMatchHistoManger::Fill_ToF(int charge, int pid, int cent, float pt, float eta, float phi)
{
  h_mTracks_ToF[charge][pid][cent]->Fill(pt,eta,phi);
  h_mTracks_ToF[charge][pid][9]->Fill(pt,eta,phi); // miniBias
  int i_eta = findEta(eta);
  int i_phi = findPhi(phi);
  string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
  h_mCounts_ToF[HistName_ToF]->Fill(pt);

  string HistName_ToF_miniBias = Form("h_mCounts_ToF_%s%s_Cent_9_Eta_%d_Phi_%d",tof::mPID_ToF[pid].c_str(),tof::mCharge[charge].c_str(),i_eta,i_phi);
  h_mCounts_ToF[HistName_ToF_miniBias]->Fill(pt);
}

void StToFMatchHistoManger::WriteHist()
{
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	h_mTracks_TPC[i_charge][i_pid][i_cent]->Write();
	h_mTracks_ToF[i_charge][i_pid][i_cent]->Write();

	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_TPC[HistName_TPC]->Write();
	    string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_ToF[HistName_ToF]->Write();
	  }
	}
      }
    }
  }
  h_FrameEta_ToF->Write();
  h_FramePhi_ToF->Write();
}
//-------------------------------------------------------------------------------------------
