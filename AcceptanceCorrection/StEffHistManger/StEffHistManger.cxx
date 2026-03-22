#include "StRoot/StEffMcPhi/StEffHistManger.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StMessMgr.h"
#include "StRoot/Utility/type.h"

ClassImp(StEffHistManger)
//
StEffHistManger::StEffHistManger(int energy, int pid, int mode, int startpt, int stoppt)
{
  mEnergy = energy;
  mMode = mode; 
  mStartPt = startpt;
  mStopPt = stoppt;

  //if( pid == 0 )
  //{
  //  mpt_first = vmsa::pt_rebin_first[mEnergy];
  //  mpt_last  = vmsa::pt_rebin_last[mEnergy];

  //  for(int i = mpt_first; i <= mpt_last; i++)
  //  {
  //    mpt_low[i] = vmsa::pt_low[mEnergy][i];
  //    mpt_up[i]  = vmsa::pt_up[mEnergy][i];
  //  }
  //}

  if( pid == 0 && mode == 0 )
  {
    mpt_first = mStartPt; //vmsa::pt_rebin_first_cent[mEnergy];
    mpt_last  = mStopPt;  //vmsa::pt_rebin_last_cent[mEnergy];

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_low_cent[mEnergy][i];
      mpt_up[i]  = vmsa::pt_up_cent[mEnergy][i];
    }
  }

  if( pid == 0 && mode == 1 )
  {
    mpt_first = mStartPt; //vmsa::pt_rebin_first_y[mEnergy];
    mpt_last  = mStopPt;  //vmsa::pt_rebin_last_y[mEnergy];

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_low_y[mEnergy][i];
      mpt_up[i]  = vmsa::pt_up_y[mEnergy][i];
    }
  }

  if( pid == 2 )
  {
    mpt_first = vmsa::pt_rebin_firstKS[mEnergy];
    mpt_last  = vmsa::pt_rebin_lastKS[mEnergy];

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_lowKS[mEnergy][i];
      mpt_up[i]  = vmsa::pt_upKS[mEnergy][i];
    }
  }
}

StEffHistManger::~StEffHistManger()
{
  /* */
}

void StEffHistManger::InitHist()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    //std::string HistName = Form("h_mMcTracks_%d",i_cent);
    //h_mMcTracks[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::y_total,-1.0*vmsa::mEtaMax,vmsa::mEtaMax);
    //h_mMcTracks[i_cent]->Sumw2();
    //HistName = Form("h_mRcTracks_%d",i_cent);
    //h_mRcTracks[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::y_total,-1.0*vmsa::mEtaMax,vmsa::mEtaMax);
    //h_mRcTracks[i_cent]->Sumw2();
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
    {
      //cos(theta*)
      std::string HistName = Form("h_mMcEffCosS_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCosSPt[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mMcEffCosSPt[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffCosS_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCosSPt[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mRcEffCosSPt[i_cent][i_pt]->Sumw2();
      //cos(theta)
      HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCosPt[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mMcEffCosPt[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCosPt[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mRcEffCosPt[i_cent][i_pt]->Sumw2();
      for(int i_y = 0; i_y < vmsa::y_total; i_y++)
      { 
        //cos(theta*)
        HistName = Form("h_mMcEffCosS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mMcEffCosSY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mMcEffCosSY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mRcEffCosS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mRcEffCosSY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mRcEffCosSY[i_cent][i_pt][i_y]->Sumw2();
        //cos(theta)
        HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mMcEffCosY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mMcEffCosY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mRcEffCosY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mRcEffCosY[i_cent][i_pt][i_y]->Sumw2();
      }
    }
  }
  //flag_eff = 0;
  //flag_eff_PtEtaPhi = 0;
  flag_eff_Cos = 0;
}

void StEffHistManger::FillHistMc(int cent, float pt, float y, float costs, float cost)
{
  //h_mMcTracks[cent]->Fill(pt,y);
  //if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  //{//20-60%
  //  h_mMcTracks[9]->Fill(pt,y,vmsa::weight[cent]);
  //}

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
    {
      if(mMode == 0) 
      {
        h_mMcEffCosSPt[cent][i_pt]->Fill(costs);
        h_mMcEffCosPt[cent][i_pt]->Fill(cost);
      }
      if(mMode == 1)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {
          if(y > vmsa::y_bin[i_y] && y < vmsa::y_bin[i_y+1])
          {
            //std::cout << "cent = " << cent << ", i_pt = " << i_pt << ", i_y = " << i_y << std::endl;
            //std::cout << "Does first plot exist?" << std::endl;
            h_mMcEffCosSY[cent][i_pt][i_y]->Fill(costs);
            //std::cout << "Yes" << std::endl;
            //std::cout << "Does second plot exist?" << std::endl;
            h_mMcEffCosY[cent][i_pt][i_y]->Fill(cost);
            //std::cout << "Yes" << std::endl;
          }
        }
      }
    }
  }
}

void StEffHistManger::FillHistRc(int cent, float pt, float y, float costs, float cost)
{
  //h_mRcTracks[cent]->Fill(pt,y);
  //if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  //{//20-60%
  //  h_mRcTracks[9]->Fill(pt,y,vmsa::weight[cent]);
  //}

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
    {
      if(mMode == 0) 
      {
        h_mRcEffCosSPt[cent][i_pt]->Fill(costs);
        h_mRcEffCosPt[cent][i_pt]->Fill(cost);
      }
      if(mMode == 1)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {
          if(y > vmsa::y_bin[i_y] && y < vmsa::y_bin[i_y+1])
          {
            h_mRcEffCosSY[cent][i_pt][i_y]->Fill(costs);
            h_mRcEffCosY[cent][i_pt][i_y]->Fill(cost);
          }
        }
      }
    }
  }
}

//float StEffHistManger::AngleShift(float phi)
//{
//  double const Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
//  double const Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};
//
//  float phi_shift = -999.0;
//  for(int psi_bin = 0; psi_bin < 3; ++psi_bin)
//  {
//    if(phi >= Psi2_low[psi_bin] && phi < Psi2_up[psi_bin])
//    {
//      phi_shift = phi - (psi_bin-1)*2.0*TMath::Pi()/2.0;
//    }
//  }
//
//  return phi_shift;
//}

TH1D* StEffHistManger::CalEffError(TH1D *h_Mc, TH1D *h_Rc, std::string HistName)
{
  TH1D* h_ratio = (TH1D*)h_Rc->Clone();
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

//void StEffHistManger::CalEfficiency()
//{
//  for(int i_cent = 0; i_cent < 10; ++i_cent)
//  {
//    std::string HistName;
//    HistName = Form("h_mMcEffPt_Cent_%d",i_cent);
//    h_mMcEffPt[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
//    HistName = Form("h_mRcEffPt_Cent_%d",i_cent);
//    h_mRcEffPt[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
//    HistName = Form("h_mEffPt_Cent_%d",i_cent);
//    h_mEffPt[i_cent] = CalEffError(h_mMcEffPt[i_cent],h_mRcEffPt[i_cent],HistName.c_str());
//
//    HistName = Form("h_mMcEffEta_Cent_%d",i_cent);
//    h_mMcEffEta[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
//    HistName = Form("h_mRcEffEta_Cent_%d",i_cent);
//    h_mRcEffEta[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
//    HistName = Form("h_mEffEta_Cent_%d",i_cent);
//    h_mEffEta[i_cent] = CalEffError(h_mMcEffEta[i_cent],h_mRcEffEta[i_cent],HistName.c_str());
//
//    HistName = Form("h_mMcEffPhi_Cent_%d",i_cent);
//    h_mMcEffPhi[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
//    HistName = Form("h_mRcEffPhi_Cent_%d",i_cent);
//    h_mRcEffPhi[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
//    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
//    h_mEffPhi[i_cent] = CalEffError(h_mMcEffPhi[i_cent],h_mRcEffPhi[i_cent],HistName.c_str());
//  }
//  flag_eff = 1;
//}

//void StEffHistManger::CalEffPtEtaPhi()
//{
//  for(int i_cent = 0; i_cent < 10; ++i_cent)
//  {
//    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
//    {
//      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
//      {
//	std::string HistNameMc = Form("h_mMcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
//	h_mMcEffPEP[HistNameMc] = (TH1D*)h_mMcTracks[i_cent]->ProjectionX(HistNameMc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);
//
//	std::string HistNameRc = Form("h_mRcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
//	h_mRcEffPEP[HistNameRc] = (TH1D*)h_mRcTracks[i_cent]->ProjectionX(HistNameRc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);
//
//	std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
//	h_mEffPEP[HistNameEff] = CalEffError(h_mMcEffPEP[HistNameMc],h_mRcEffPEP[HistNameRc],HistNameEff.c_str());
//      }
//    }
//  }
//  flag_eff_PtEtaPhi = 1;
//}

void StEffHistManger::CalEffCosThetaStar()
{
  int cent_bin = 0;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
    {
      if(mMode == 0)
      {
        std::string HistName = Form("h_mEffCosS_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosSPt[i_cent][i_pt] = (TH1D*)CalEffError(h_mMcEffCosSPt[i_cent][i_pt],h_mRcEffCosSPt[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosPt[i_cent][i_pt] = (TH1D*)CalEffError(h_mMcEffCosPt[i_cent][i_pt],h_mRcEffCosPt[i_cent][i_pt],HistName.c_str());
      }
      if(mMode == 1)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {      
          std::string HistNameStar = Form("h_mEffCosS_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);
          std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);

          if(i_cent > vmsa::cent_rebin[cent_bin] && i_cent < vmsa::cent_rebin[cent_bin])
          {
            h_mMcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffCosSY[i_cent][i_pt][i_y],1.0);
            h_mRcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffCosSY[i_cent][i_pt][i_y],1.0);
          }
          if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) 
          {
            h_mEffCosSY[cent_bin][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistNameStar.c_str());
            h_mEffCosY[cent_bin][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistName.c_str());
          }
        }
      }
    }
    if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) cent_bin++;
  }
  flag_eff_Cos = 1;
}

void StEffHistManger::WriteHist()
{
  if(flag_eff_Cos > 0.5)
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
      {
        if(mMode == 0)
        {
          h_mMcEffCosSPt[i_cent][i_pt]->Write();
          h_mRcEffCosSPt[i_cent][i_pt]->Write();
          h_mEffCosSPt[i_cent][i_pt]->Write();
          h_mMcEffCosPt[i_cent][i_pt]->Write();
          h_mRcEffCosPt[i_cent][i_pt]->Write();
          h_mEffCosPt[i_cent][i_pt]->Write();
        }
        if(mMode == 1)
        {
          for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
          {
            h_mMcEffCosSY[i_cent][i_pt][i_y]->Write();
            h_mRcEffCosSY[i_cent][i_pt][i_y]->Write();
            h_mEffCosSY[i_cent][i_pt][i_y]->Write();
            h_mMcEffCosY[i_cent][i_pt][i_y]->Write();
            h_mRcEffCosY[i_cent][i_pt][i_y]->Write();
            h_mEffCosY[i_cent][i_pt][i_y]->Write();
          }
        }
      }
    }
  }
}
