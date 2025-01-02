#include "StEffHistManger.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include <iostream>
#include "StRoot/Utility/StSpinAlignmentCons.h"
// #include <string>

ClassImp(StEffHistManger)
//
StEffHistManger::StEffHistManger()
{
  /* */
}

StEffHistManger::~StEffHistManger()
{
  /* */
}

void StEffHistManger::InitHist()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName = Form("h_mMcTracks_%d",i_cent);
    h_mMcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,vmsa::ptMin,vmsa::ptMax,vmsa::BinEta,-vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,-TMath::Pi(),TMath::Pi());
    //h_mMcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,-1.0,1.0,vmsa::BinPhi,-0.5*TMath::Pi(),0.5*TMath::Pi());
    h_mMcTracks[i_cent]->Sumw2();
    HistName = Form("h_mRcTracks_%d",i_cent);
    //h_mRcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,-1.0,1.0,vmsa::BinPhi,-0.5*TMath::Pi(),0.5*TMath::Pi());
    h_mRcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,vmsa::ptMin,vmsa::ptMax,vmsa::BinEta,-vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,-TMath::Pi(),TMath::Pi());
    h_mRcTracks[i_cent]->Sumw2();
    HistName = Form("h_mPtGl_%d",i_cent);
    h_mPtGl[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinPt*10,0.0,24.0);
    HistName = Form("h_mPtPr_%d",i_cent);
    h_mPtPr[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinPt*10,0.0,24.0);
  }
  h_FrameEta = new TH1D("h_FrameEta","h_FrameEta",vmsa::BinEta,vmsa::mEtaMax,vmsa::mEtaMax);
  h_FramePhi = new TH1D("h_FramePhi","h_FramePhi",vmsa::BinPhi,-TMath::Pi(),TMath::Pi());
  //h_FramePhi = new TH1D("h_FramePhi","h_FramePhi",vmsa::BinPhi,-0.5*TMath::Pi(),0.5*TMath::Pi());
  flag_eff = 0;
  flag_eff_PtEtaPhi = 0;
}

float StEffHistManger::AngleShift(float phi)
{
  double const Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
  double const Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};

  float phi_shift = -999.0;
  for(int psi_bin = 0; psi_bin < 3; ++psi_bin)
  {
    if(phi >= Psi2_low[psi_bin] && phi < Psi2_up[psi_bin])
    {
      phi_shift = phi - (psi_bin-1)*2.0*TMath::Pi()/2.0;
    }
  }

  return phi_shift;
}

void StEffHistManger::FillHistMc(int cent, float pt, float eta, float phi, float psiWest, float psiEast)
{
  h_mMcTracks[cent]->Fill(pt,eta,phi);
  if(cent >= 2 && cent <= 5) h_mMcTracks[9]->Fill(pt,eta,phi);
  //if(eta > -1.0*vmsa::mEtaMax && eta < 0.0)
  //
  //{ // east eta cut
  //  float phi_shift = AngleShift(phi-psiWest);
  //  h_mMcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mMcTracks[9]->Fill(pt,eta,phi_shift);
  //}
  //if(eta >= 0.0 && eta < vmsa::mEtaMax)
  //{ // west eta cut
  //  float phi_shift = AngleShift(phi-psiEast);
  //  h_mMcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mMcTracks[9]->Fill(pt,eta,phi_shift);
  //}
}

void StEffHistManger::FillHistRc(int cent, float pt, float eta, float phi, float psiWest, float psiEast)
{
  h_mRcTracks[cent]->Fill(pt,eta,phi);
  if(cent >= 2 && cent <= 5) h_mRcTracks[9]->Fill(pt,eta,phi);
  //if(eta > -1.0*vmsa::mEtaMax && eta < 0.0)
  //{ // east eta cut
  //  float phi_shift = AngleShift(phi-psiWest);
  //  h_mRcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mRcTracks[9]->Fill(pt,eta,phi_shift);
  //}
  //if(eta >= 0.0 && eta < vmsa::mEtaMax)
  //{ // west eta cut
  //  float phi_shift = AngleShift(phi-psiEast);
  //  h_mRcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mRcTracks[9]->Fill(pt,eta,phi_shift);
  //}
}

void StEffHistManger::FillHistPt(int cent, float McPt, float gRcPt, float pRcPt)
{
  h_mPtGl[cent]->Fill(McPt,gRcPt-McPt);
  if(cent >= 2 && cent <= 5) h_mPtGl[9]->Fill(McPt,gRcPt-McPt);
  h_mPtPr[cent]->Fill(McPt,pRcPt-McPt);
  if(cent >= 2 && cent <= 5) h_mPtPr[9]->Fill(McPt,pRcPt-McPt);
}

TH1D* StEffHistManger::CalEffError(TH1D *h_Mc, TH1D *h_Rc, std::string HistName)
{
  TH1D* h_ratio = (TH1D*)h_Rc->Clone();
  //h_ratio->Divide(h_Mc);
  h_ratio->Divide(h_Rc,h_Mc,1,1,"B");
  //for(int i_bin = 1; i_bin < h_ratio->GetNbinsX()+1; ++i_bin)
  //{
  //  double n = h_Mc->GetBinContent(i_bin);
  //  double k = h_Rc->GetBinContent(i_bin);
  //  double variance = (k+1.0)*(k+2.0)/((n+2.0)*(n+3.0))-(k+1.0)*(k+1.0)/((n+2.0)*(n+2.0));
  //  double sigma = TMath::Sqrt(variance);
  //  if(n > 0.0 && k > 0.0) h_ratio->SetBinError(i_bin,sigma);
  //}
  h_ratio->SetName(HistName.c_str());

  return h_ratio;
}

void StEffHistManger::CalEfficiency()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName;

    HistName = Form("h_mMcEffPt_Cent_%d",i_cent);
    h_mMcEffPt[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPt_Cent_%d",i_cent);
    h_mRcEffPt[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffPt[i_cent] = CalEffError(h_mMcEffPt[i_cent],h_mRcEffPt[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffEta_Cent_%d",i_cent);
    h_mMcEffEta[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffEta_Cent_%d",i_cent);
    h_mRcEffEta[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEffEta[i_cent] = CalEffError(h_mMcEffEta[i_cent],h_mRcEffEta[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffPhi_Cent_%d",i_cent);
    h_mMcEffPhi[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPhi_Cent_%d",i_cent);
    h_mRcEffPhi[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffPhi[i_cent] = CalEffError(h_mMcEffPhi[i_cent],h_mRcEffPhi[i_cent],HistName.c_str());
  }
  flag_eff = 1;
}

void StEffHistManger::CalEffPtEtaPhi()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	std::string HistNameMc = Form("h_mMcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mMcEffPEP[HistNameMc] = (TH1D*)h_mMcTracks[i_cent]->ProjectionX(HistNameMc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);

	std::string HistNameRc = Form("h_mRcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mRcEffPEP[HistNameRc] = (TH1D*)h_mRcTracks[i_cent]->ProjectionX(HistNameRc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);

	std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEffPEP[HistNameEff] = CalEffError(h_mMcEffPEP[HistNameMc],h_mRcEffPEP[HistNameRc],HistNameEff.c_str());
      }
    }
  }
  flag_eff_PtEtaPhi = 1;
}

void StEffHistManger::WriteHist()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    h_mMcTracks[i_cent]->Write();
    h_mRcTracks[i_cent]->Write();
    h_mPtGl[i_cent]->Write();
    h_mPtPr[i_cent]->Write();
    //if(flag_eff < 1) continue;
    //h_mEffPt[i_cent]->Write();
    //h_mEffEta[i_cent]->Write();
    //h_mEffPhi[i_cent]->Write();

    //if(flag_eff_PtEtaPhi < 1) continue;
    //for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    //{
    //  for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
    //  {
    //    std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
    //    h_mEffPEP[HistNameEff]->Write();
    //  }
    //}
  }
  h_FrameEta->Write();
  h_FramePhi->Write();
}
