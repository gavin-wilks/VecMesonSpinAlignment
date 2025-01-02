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

  if( pid == 0 && mode == 2 )
  {
    mpt_first = mStartPt; //vmsa::pt_rebin_first_y[mEnergy];
    mpt_last  = mStopPt;  //vmsa::pt_rebin_last_y[mEnergy];

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_low[mEnergy][i];
      mpt_up[i]  = vmsa::pt_up[mEnergy][i];
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
      HistName = Form("p_mMcCos2BetaP_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mMcCos2BetaPPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mRcCos2BetaP_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mRcCos2BetaPPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mMcCos4BetaP_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mMcCos4BetaPPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mRcCos4BetaP_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mRcCos4BetaPPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mMcCos2BetaPCos4BetaP_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mMcCos2BetaPCos4BetaPPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mRcCos2BetaPCos4BetaP_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mRcCos2BetaPCos4BetaPPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      //cos(theta*) RP
      HistName = Form("h_mMcEffCosS_RP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCosSPt_RP[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mMcEffCosSPt_RP[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffCosS_RP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCosSPt_RP[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mRcEffCosSPt_RP[i_cent][i_pt]->Sumw2();
      HistName = Form("p_mMcCos2Beta_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mMcCos2BetaPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mRcCos2Beta_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mRcCos2BetaPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mMcCos4Beta_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mMcCos4BetaPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mRcCos4Beta_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mRcCos4BetaPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mMcCos2BetaCos4Beta_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mMcCos2BetaCos4BetaPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      HistName = Form("p_mRcCos2BetaCos4Beta_Cent_%d_Pt_%d",i_cent,i_pt);
      p_mRcCos2BetaCos4BetaPt[i_cent][i_pt] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
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
        HistName = Form("p_mMcCos2BetaP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mMcCos2BetaPY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mRcCos2BetaP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mRcCos2BetaPY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mMcCos4BetaP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mMcCos4BetaPY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mRcCos4BetaP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mRcCos4BetaPY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mMcCos2BetaPCos4BetaP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mMcCos2BetaPCos4BetaPY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mRcCos2BetaPCos4BetaP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mRcCos2BetaPCos4BetaPY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        //cos(theta*) w.r.t. RP
        HistName = Form("h_mMcEffCosS_RP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mMcEffCosSY_RP[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mMcEffCosSY_RP[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mRcEffCosS_RP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mRcEffCosSY_RP[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mRcEffCosSY_RP[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("p_mMcCos2Beta_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mMcCos2BetaY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mRcCos2Beta_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mRcCos2BetaY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mMcCos4Beta_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mMcCos4BetaY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mRcCos4Beta_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mRcCos4BetaY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mMcCos2BetaCos4Beta_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mMcCos2BetaCos4BetaY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        HistName = Form("p_mRcCos2BetaCos4Beta_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        p_mRcCos2BetaCos4BetaY[i_cent][i_pt][i_y] = new TProfile(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
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

void StEffHistManger::FillHistMc(int cent, float pt, float y, float cost, float costs, float costsRP, float betaP, float beta)
{
  //h_mMcTracks[cent]->Fill(pt,y);
  //if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  //{//20-60%
  //  h_mMcTracks[9]->Fill(pt,y,vmsa::weight[cent]);
  //}
  double cos2b        = TMath::Cos(2.0*beta);
  double cos2bP       = TMath::Cos(2.0*betaP);
  double cos4b        = TMath::Cos(4.0*beta);
  double cos4bP       = TMath::Cos(4.0*betaP);
  double cos2bcos4b   = cos2b*cos4b;
  double cos2bPcos4bP = cos2bP*cos4bP;

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
    {
      if(mMode == 0) 
      {
        h_mMcEffCosSPt[cent][i_pt]->Fill(costs);
        h_mMcEffCosSPt_RP[cent][i_pt]->Fill(costsRP);
        h_mMcEffCosPt[cent][i_pt]->Fill(cost);

        p_mMcCos2BetaPPt[cent][i_pt]->Fill(costs,cos2bP);
        p_mMcCos2BetaPt[cent][i_pt]->Fill(costsRP,cos2b);
        p_mMcCos4BetaPPt[cent][i_pt]->Fill(costs,cos4bP);
        p_mMcCos4BetaPt[cent][i_pt]->Fill(costsRP,cos4b);
        p_mMcCos2BetaPCos4BetaPPt[cent][i_pt]->Fill(costs,cos2bPcos4bP);
        p_mMcCos2BetaCos4BetaPt[cent][i_pt]->Fill(costsRP,cos2bcos4b);
      }
      if(mMode == 1)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {
          if(y > vmsa::y_bin[i_y] && y < vmsa::y_bin[i_y+1])
          {
            h_mMcEffCosSY[cent][i_pt][i_y]->Fill(costs);
            h_mMcEffCosSY_RP[cent][i_pt][i_y]->Fill(costsRP);
            h_mMcEffCosY[cent][i_pt][i_y]->Fill(cost);

            p_mMcCos2BetaPY[cent][i_pt][i_y]->Fill(costs,cos2bP);
            p_mMcCos2BetaY[cent][i_pt][i_y]->Fill(costsRP,cos2b);
            p_mMcCos4BetaPY[cent][i_pt][i_y]->Fill(costs,cos4bP);
            p_mMcCos4BetaY[cent][i_pt][i_y]->Fill(costsRP,cos4b);
            p_mMcCos2BetaPCos4BetaPY[cent][i_pt][i_y]->Fill(costs,cos2bPcos4bP);
            p_mMcCos2BetaCos4BetaY[cent][i_pt][i_y]->Fill(costsRP,cos2bcos4b);
          }
        }
      }
      if(mMode == 2)
      {
        for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
        {
          if(y > vmsa::y_pt_rebin[i_y] && y < vmsa::y_pt_rebin[i_y+1])
          {
            //std::cout << "cent = " << cent << ", i_pt = " << i_pt << ", i_y = " << i_y << std::endl;
            //std::cout << "Does first plot exist?" << std::endl;
            h_mMcEffCosSY[cent][i_pt][i_y]->Fill(costs);
            h_mMcEffCosSY_RP[cent][i_pt][i_y]->Fill(costsRP);
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

void StEffHistManger::FillHistRc(int cent, float pt, float y, float cost, float costs, float costsRP, float betaP, float beta)
{
  //h_mRcTracks[cent]->Fill(pt,y);
  //if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  //{//20-60%
  //  h_mRcTracks[9]->Fill(pt,y,vmsa::weight[cent]);
  //}
  double cos2b        = TMath::Cos(2.0*beta);
  double cos2bP       = TMath::Cos(2.0*betaP);
  double cos4b        = TMath::Cos(4.0*beta);
  double cos4bP       = TMath::Cos(4.0*betaP);
  double cos2bcos4b   = cos2b*cos4b;
  double cos2bPcos4bP = cos2bP*cos4bP;

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
    {
      if(mMode == 0) 
      {
        h_mRcEffCosSPt[cent][i_pt]->Fill(costs);
        h_mRcEffCosSPt_RP[cent][i_pt]->Fill(costsRP);
        h_mRcEffCosPt[cent][i_pt]->Fill(cost);
        
        p_mRcCos2BetaPPt[cent][i_pt]->Fill(costs,cos2bP);
        p_mRcCos2BetaPt[cent][i_pt]->Fill(costsRP,cos2b);
        p_mRcCos4BetaPPt[cent][i_pt]->Fill(costs,cos4bP);
        p_mRcCos4BetaPt[cent][i_pt]->Fill(costsRP,cos4b);
        p_mRcCos2BetaPCos4BetaPPt[cent][i_pt]->Fill(costs,cos2bPcos4bP);
        p_mRcCos2BetaCos4BetaPt[cent][i_pt]->Fill(costsRP,cos2bcos4b);
      }
      if(mMode == 1)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {
          if(y > vmsa::y_bin[i_y] && y < vmsa::y_bin[i_y+1])
          {
            h_mRcEffCosSY[cent][i_pt][i_y]->Fill(costs);
            h_mRcEffCosSY_RP[cent][i_pt][i_y]->Fill(costsRP);
            h_mRcEffCosY[cent][i_pt][i_y]->Fill(cost);

            p_mRcCos2BetaPY[cent][i_pt][i_y]->Fill(costs,cos2bP);
            p_mRcCos2BetaY[cent][i_pt][i_y]->Fill(costsRP,cos2b);
            p_mRcCos4BetaPY[cent][i_pt][i_y]->Fill(costs,cos4bP);
            p_mRcCos4BetaY[cent][i_pt][i_y]->Fill(costsRP,cos4b);
            p_mRcCos2BetaPCos4BetaPY[cent][i_pt][i_y]->Fill(costs,cos2bPcos4bP);
            p_mRcCos2BetaCos4BetaY[cent][i_pt][i_y]->Fill(costsRP,cos2bcos4b);
          }
        }
      }
      if(mMode == 2)
      {
        for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
        {
          if(y > vmsa::y_pt_rebin[i_y] && y < vmsa::y_pt_rebin[i_y+1])
          {
            //std::cout << "cent = " << cent << ", i_pt = " << i_pt << ", i_y = " << i_y << std::endl;
            //std::cout << "Does first plot exist?" << std::endl;
            h_mRcEffCosSY[cent][i_pt][i_y]->Fill(costs);
            h_mRcEffCosSY_RP[cent][i_pt][i_y]->Fill(costsRP);
            //std::cout << "Yes" << std::endl;
            //std::cout << "Does second plot exist?" << std::endl;
            h_mRcEffCosY[cent][i_pt][i_y]->Fill(cost);
            //std::cout << "Yes" << std::endl;
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
  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
  {
    if(mMode == 2)
    {
      for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
      {      
        std::string HistNameStar = Form("h_mEffCosS_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
        std::string HistNameStarRP = Form("h_mEffCosS_RP_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
        std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);

        h_mEffCosSY[9][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosSY[9][i_pt][i_y],h_mRcEffCosSY[9][i_pt][i_y],HistNameStar.c_str());
        h_mEffCosSY_RP[9][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosSY_RP[9][i_pt][i_y],h_mRcEffCosSY_RP[9][i_pt][i_y],HistNameStarRP.c_str());
        h_mEffCosY[9][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosY[9][i_pt][i_y],h_mRcEffCosY[9][i_pt][i_y],HistName.c_str());
      }
    }
  }



  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
    {
      if(mMode == 0)
      {
        std::string HistName = Form("h_mEffCosS_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosSPt[i_cent][i_pt] = (TH1D*)CalEffError(h_mMcEffCosSPt[i_cent][i_pt],h_mRcEffCosSPt[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCosS_RP_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosSPt_RP[i_cent][i_pt] = (TH1D*)CalEffError(h_mMcEffCosSPt_RP[i_cent][i_pt],h_mRcEffCosSPt_RP[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosPt[i_cent][i_pt] = (TH1D*)CalEffError(h_mMcEffCosPt[i_cent][i_pt],h_mRcEffCosPt[i_cent][i_pt],HistName.c_str());

        std::string GraphNameStar = Form("g_mEffCosS_Cent_%d_Pt_%d",i_cent,i_pt);
        g_mEffCosSPt[i_cent][i_pt] = new TGraphAsymmErrors();
        g_mEffCosSPt[i_cent][i_pt]->BayesDivide(h_mRcEffCosSPt[i_cent][i_pt],h_mMcEffCosSPt[i_cent][i_pt]);
        g_mEffCosSPt[i_cent][i_pt]->SetName(GraphNameStar.c_str());

        std::string GraphName = Form("g_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
        g_mEffCosPt[i_cent][i_pt] = new TGraphAsymmErrors();
        g_mEffCosPt[i_cent][i_pt]->BayesDivide(h_mRcEffCosPt[i_cent][i_pt],h_mMcEffCosPt[i_cent][i_pt]);
        g_mEffCosPt[i_cent][i_pt]->SetName(GraphName.c_str());
      }
      if(mMode == 1)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {      
          std::string HistNameStar = Form("h_mEffCosS_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);
          std::string HistNameStarRP = Form("h_mEffCosS_RP_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);
          std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);
          std::string GraphNameStar = Form("g_mEffCosS_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);
          std::string GraphName = Form("g_mEffCos_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);

          if(i_cent > vmsa::cent_rebin[cent_bin] && i_cent < vmsa::cent_rebin[cent_bin+1])
          {
            h_mMcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffCosSY[i_cent][i_pt][i_y],1.0);
            h_mMcEffCosSY_RP[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffCosSY_RP[i_cent][i_pt][i_y],1.0);
            h_mRcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffCosSY[i_cent][i_pt][i_y],1.0);
            h_mRcEffCosSY_RP[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffCosSY_RP[i_cent][i_pt][i_y],1.0);
            h_mMcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffCosY[i_cent][i_pt][i_y],1.0);
            h_mRcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffCosY[i_cent][i_pt][i_y],1.0);
          }
          if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) 
          {
            h_mEffCosSY[cent_bin][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistNameStar.c_str());
            h_mEffCosSY_RP[cent_bin][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosSY_RP[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosSY_RP[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistNameStarRP.c_str());
            h_mEffCosY[cent_bin][i_pt][i_y] = (TH1D*)CalEffError(h_mMcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistName.c_str());
            g_mEffCosSY[cent_bin][i_pt][i_y] = new TGraphAsymmErrors();
            g_mEffCosSY[cent_bin][i_pt][i_y]->BayesDivide(h_mRcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mMcEffCosSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]);
            g_mEffCosSY[cent_bin][i_pt][i_y]->SetName(GraphNameStar.c_str());
            g_mEffCosY[cent_bin][i_pt][i_y] = new TGraphAsymmErrors();
            g_mEffCosY[cent_bin][i_pt][i_y]->BayesDivide(h_mRcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mMcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]);
            g_mEffCosY[cent_bin][i_pt][i_y]->SetName(GraphName.c_str());
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
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
    {
      if(mMode == 2)
      {
        for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
        {      
          std::string HistName = Form("h_mMcEffCosS_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
          h_mMcEffCosSY[9][i_pt][i_y]->SetName(HistName.c_str());
          h_mMcEffCosSY[9][i_pt][i_y]->Write();
          HistName = Form("h_mMcEffCosS_RP_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
          h_mMcEffCosSY_RP[9][i_pt][i_y]->SetName(HistName.c_str());
          h_mMcEffCosSY_RP[9][i_pt][i_y]->Write();
          HistName = Form("h_mRcEffCosS_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
          h_mRcEffCosSY[9][i_pt][i_y]->SetName(HistName.c_str());
          h_mRcEffCosSY[9][i_pt][i_y]->Write();
          HistName = Form("h_mRcEffCosS_RP_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
          h_mRcEffCosSY_RP[9][i_pt][i_y]->SetName(HistName.c_str());
          h_mRcEffCosSY_RP[9][i_pt][i_y]->Write();
          h_mEffCosSY_RP[9][i_pt][i_y]->Write();
          h_mEffCosSY[9][i_pt][i_y]->Write();
          HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
          h_mMcEffCosY[9][i_pt][i_y]->SetName(HistName.c_str());
          h_mMcEffCosY[9][i_pt][i_y]->Write();
          HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d_Y_%d",9,i_pt,i_y);
          h_mRcEffCosY[9][i_pt][i_y]->SetName(HistName.c_str());
          h_mRcEffCosY[9][i_pt][i_y]->Write();
          h_mEffCosY[9][i_pt][i_y]->Write(); 
        }
      }
    }
  
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
      {
        if(mMode == 0)
        {
          h_mMcEffCosSPt[i_cent][i_pt]->Write();
          h_mRcEffCosSPt[i_cent][i_pt]->Write();
          h_mEffCosSPt[i_cent][i_pt]->Write();
          h_mMcEffCosSPt_RP[i_cent][i_pt]->Write();
          h_mRcEffCosSPt_RP[i_cent][i_pt]->Write();
          h_mEffCosSPt_RP[i_cent][i_pt]->Write();
          h_mMcEffCosPt[i_cent][i_pt]->Write();
          h_mRcEffCosPt[i_cent][i_pt]->Write();
          h_mEffCosPt[i_cent][i_pt]->Write();
  
          p_mMcCos2BetaPPt[i_cent][i_pt]->Write();
          p_mMcCos2BetaPt[i_cent][i_pt]->Write();
          p_mMcCos4BetaPPt[i_cent][i_pt]->Write();
          p_mMcCos4BetaPt[i_cent][i_pt]->Write();
          p_mMcCos2BetaPCos4BetaPPt[i_cent][i_pt]->Write();
          p_mMcCos2BetaCos4BetaPt[i_cent][i_pt]->Write();
          p_mRcCos2BetaPPt[i_cent][i_pt]->Write();
          p_mRcCos2BetaPt[i_cent][i_pt]->Write();
          p_mRcCos4BetaPPt[i_cent][i_pt]->Write();
          p_mRcCos4BetaPt[i_cent][i_pt]->Write();
          p_mRcCos2BetaPCos4BetaPPt[i_cent][i_pt]->Write();
          p_mRcCos2BetaCos4BetaPt[i_cent][i_pt]->Write();

          g_mEffCosSPt[i_cent][i_pt]->Write();
          g_mEffCosPt[i_cent][i_pt]->Write();
        }
        if(mMode == 1)
        { 
          if( i_cent >= vmsa::cent_rebin_total ) continue;
          for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
          {
            std::string HistName = Form("h_mMcEffCosS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
            h_mMcEffCosSY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->SetName(HistName.c_str());
            h_mMcEffCosSY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            HistName = Form("h_mMcEffCosS_RP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
            h_mMcEffCosSY_RP[vmsa::cent_rebin[i_cent]][i_pt][i_y]->SetName(HistName.c_str());
            h_mMcEffCosSY_RP[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            HistName = Form("h_mRcEffCosS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
            h_mRcEffCosSY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->SetName(HistName.c_str());
            h_mRcEffCosSY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            HistName = Form("h_mRcEffCosS_RP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
            h_mRcEffCosSY_RP[vmsa::cent_rebin[i_cent]][i_pt][i_y]->SetName(HistName.c_str());
            h_mRcEffCosSY_RP[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            h_mEffCosSY_RP[i_cent][i_pt][i_y]->Write();
            h_mEffCosSY[i_cent][i_pt][i_y]->Write();
            HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
            h_mMcEffCosY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->SetName(HistName.c_str());
            h_mMcEffCosY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
            h_mRcEffCosY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->SetName(HistName.c_str());
            h_mRcEffCosY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            h_mEffCosY[i_cent][i_pt][i_y]->Write();

            p_mMcCos2BetaPY[i_cent][i_pt][i_y]->Write();
            p_mMcCos2BetaY[i_cent][i_pt][i_y]->Write();
            p_mMcCos4BetaPY[i_cent][i_pt][i_y]->Write();
            p_mMcCos4BetaY[i_cent][i_pt][i_y]->Write();
            p_mMcCos2BetaPCos4BetaPY[i_cent][i_pt][i_y]->Write();
            p_mMcCos2BetaCos4BetaY[i_cent][i_pt][i_y]->Write();
            p_mRcCos2BetaPY[i_cent][i_pt][i_y]->Write();
            p_mRcCos2BetaY[i_cent][i_pt][i_y]->Write();
            p_mRcCos4BetaPY[i_cent][i_pt][i_y]->Write();
            p_mRcCos4BetaY[i_cent][i_pt][i_y]->Write();
            p_mRcCos2BetaPCos4BetaPY[i_cent][i_pt][i_y]->Write();
            p_mRcCos2BetaCos4BetaY[i_cent][i_pt][i_y]->Write();
           
            g_mEffCosSY[i_cent][i_pt][i_y]->Write();
            g_mEffCosY[i_cent][i_pt][i_y]->Write();
          }
        }
      }
    }
  }
}
