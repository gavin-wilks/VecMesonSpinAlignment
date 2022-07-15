#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calFlowKStar_Poly3(int energy = 4, int pid = 2, int year = 0)
{
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptDate(0);

  string inputfile = Form("../data/file_%s_Resolution_20220330.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Res = TFile::Open(inputfile.c_str());

  double mTpcSubRes2Val[9];
  double mTpcSubRes2Err[9];
  double mTpcFullRes2Val[9];
  double mTpcFullRes2Err[9];

  // calculate sub event plane resolution
  TProfile *p_mTpcSubRes2 = (TProfile*)File_Res->Get("p_mRes2_Sub");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcSubRes2->GetBinContent(p_mTpcSubRes2->FindBin(i_cent));
    const double errRaw = p_mTpcSubRes2->GetBinError(p_mTpcSubRes2->FindBin(i_cent));
    if(resRaw > 0)
    {
      mTpcSubRes2Val[i_cent] = TMath::Sqrt(resRaw);
      mTpcSubRes2Err[i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));
    }
  }

  string InPutFile_SE = "../data/Yields_KStar_SE_19GeV_20220708.root";
  
  string InPutFile_ME = "../data/Yields_KStar_ME_19GeV_20220708.root";
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(int i_pt = vmsa::pt_startKS; i_pt < vmsa::pt_stopKS; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
    {
      for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
      {
        string KEY_SE = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
        string KEY_ME = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
           
        h_mMass_SE[KEY_SE] = (TH1F*)File_SE->Get(KEY_SE.c_str())->Clone(); 
        int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][0]);
        int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][0]);
        float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
      
        h_mMass_ME[KEY_ME] = (TH1F*)File_ME->Get(KEY_ME.c_str())->Clone();
        float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
        h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
      
        string KEY_SM = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
        h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
        h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
        cout << KEY_SE << endl;
      }
    }
  }
  cout << "Out of the first loop" << endl;
#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *cy6 = new TCanvas("cy6","cy6",10,10,800,600);
  cy6->cd();
  string KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",6,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",6,9,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",6,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy6->SaveAs("./figures/cy6_costheta3_pt6_cent9.pdf");

  TCanvas *cy7 = new TCanvas("cy7","cy7",10,10,800,600);
  cy7->cd();
  KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",7,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",7,9,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",7,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy7->SaveAs("./figures/cy7_costheta3_pt7_cent9.pdf");

  TCanvas *cy5 = new TCanvas("cy5","cy5",10,10,800,600);
  cy5->cd();
  KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",5,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",5,9,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",5,9,0,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy5->SaveAs("./figures/cy5_costheta3_pt5_cent9.pdf");
 


  // QA Plots for pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
  c_pT->Divide(5,5);
  for(int i_pt = vmsa::pt_startKS; i_pt < vmsa::pt_stopKS; i_pt++) // pt loop
  {
    c_pT->cd(i_pt+1);
    c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT->cd(i_pt+1)->SetTicks(1,1);
    c_pT->cd(i_pt+1)->SetGrid(0,0);
    string KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",i_pt,9,3,vmsa::mPID[pid].c_str());
    h_mMass_SE[KEY_SE_QA]->DrawCopy();

    string KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",i_pt,9,3,vmsa::mPID[pid].c_str());
    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,9,3,vmsa::mPID[pid].c_str());
    h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStartKS[i_pt],vmsa::ptRawStopKS[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
  {
    for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
    {
      for(int pt_binKS = vmsa::pt_rebin_firstKSv2[energy]; pt_binKS < vmsa::pt_rebin_lastKSv2[energy]; pt_binKS++) // pt loop
      {
        string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",pt_binKS,i_cent,i_theta,vmsa::mPID[pid].c_str());
        for(int i_pt = vmsa::pt_rebin_startKSv2[energy][pt_binKS]; i_pt <= vmsa::pt_rebin_stopKSv2[energy][pt_binKS]; i_pt++)
        {
          string KEY_SM = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
          if(i_pt == vmsa::pt_rebin_startKSv2[energy][pt_binKS])
          {
            h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone();
          }
          else
          {
            h_mMass[KEY]->Add(h_mMass_SM[KEY_SM],1.0);
          }
        }
      }
    }
  }

#if _PlotQA_
  // QA Plots for pT rebins
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
  {
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,9,3,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetMarkerStyle(20);
    h_mMass[KEY_QA]->SetMarkerSize(0.4);
    h_mMass[KEY_QA]->SetLineColor(1);
    h_mMass[KEY_QA]->DrawCopy("pE");

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

    // Poly + Breit Wignar fit to phi integrated InvMass
    TH1FMap h_mMass_theta;
    vecFMap ParFit_theta;

    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
      {
	string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
          //h_mMass[KEY]->Rebin(2);
	  if(i_theta == 0) h_mMass_theta[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone();
	  else h_mMass_theta[KEY_theta]->Add(h_mMass[KEY],1.0);
	}
	TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	for(int i_par = 0; i_par < 7; i_par++)
	{
	  f_bw->ReleaseParameter(i_par);
	}
	f_bw->SetParameter(0,vmsa::InvMass[pid]);
	f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.01,vmsa::InvMass[pid]+0.01);
	f_bw->SetParameter(1,vmsa::Width[pid]);
	f_bw->SetParameter(2,10000);
	//f_bw->SetParameter(3,6000);
	//f_bw->SetParameter(4,0.5);
        //f_bw->SetParameter(5,0.1);
	//f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	ParFit_theta[KEY_theta].clear();
	h_mMass_theta[KEY_theta]->Fit(f_bw,"NQR");
	for(int n_par = 0; n_par < 7; n_par++)
	{
	  ParFit_theta[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(n_par)));
	}
      }
    }

#if _PlotQA_
    // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
    TCanvas *c_linsub = new TCanvas("c_linsub","c_linsub",10,10,800,600);
    TCanvas *c_mMass_theta = new TCanvas("c_mMass_theta","c_mMass_theta",10,10,1400,1400);
    c_mMass_theta->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,9,vmsa::mPID[pid].c_str()); 
      h_mMass_theta[KEY_theta_QA]->SetMarkerColor(1);
      h_mMass_theta[KEY_theta_QA]->SetMarkerStyle(24);
      h_mMass_theta[KEY_theta_QA]->SetMarkerSize(0.8);
      h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");

      TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      for(int i_par = 0; i_par < 7; i_par++)
      {
	f_bw->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par]);
      }
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->DrawCopy("l same");

      TF1 *f_poly = new TF1("f_poly",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      f_poly->SetParameter(0,ParFit_theta[KEY_theta_QA][0]);
      f_poly->SetParameter(1,ParFit_theta[KEY_theta_QA][1]);
      f_poly->SetParameter(2,ParFit_theta[KEY_theta_QA][2]);
      f_poly->SetParameter(3,ParFit_theta[KEY_theta_QA][3]);
      f_poly->SetParameter(4,ParFit_theta[KEY_theta_QA][4]);
      f_poly->SetParameter(5,ParFit_theta[KEY_theta_QA][5]);
      f_poly->SetParameter(6,ParFit_theta[KEY_theta_QA][6]);
      f_poly->SetLineColor(4);
      f_poly->SetLineStyle(2);
      f_poly->SetLineWidth(4);
      f_poly->DrawCopy("l same");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);

      if(i_pt == 2)
      {
        c_linsub->cd();
        h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");
        f_bw->DrawCopy("l same");
        f_poly->DrawCopy("l same");

        string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
        plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
        PlotLine(0.98,1.05,0.0,0.0,1,2,2);
        c_linsub->SaveAs("./figures/c_linsub.pdf");
      }
    }
    c_mMass_theta->SaveAs("./figures/c_mMass_theta.pdf");
#endif

    // Poly+bw fits for phi differential InvMass
    //for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    //{
    //  for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
    //  {
    //    string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
    //    for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
    //    {
    //      string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
    //      TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6); //Poly+bw fits
    //      for(int i_par = 0; i_par < 6; i_par++)
    //      {
    //        f_bw->ReleaseParameter(i_par);
    //      }
    //      f_bw->FixParameter(0,ParFit_theta[KEY_theta][0]);
    //      f_bw->FixParameter(1,ParFit_theta[KEY_theta][1]);
    //      f_bw->SetParameter(2,ParFit_theta[KEY_theta][2]/10.0);
    //      f_bw->SetParameter(3,ParFit_theta[KEY_theta][3]/10.0);
    //      f_bw->SetParameter(4,ParFit_theta[KEY_theta][4]/10.0);
    //      f_bw->SetParameter(5,ParFit_theta[KEY_theta][5]/10.0);
    //      f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    //      h_mMass[KEY]->Fit(f_bw,"NQR");

    //      TF1 *f_poly = new TF1("f_poly",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
    //      f_poly->FixParameter(0,f_bw->GetParameter(3));
    //      f_poly->FixParameter(1,f_bw->GetParameter(4));
    //      f_poly->FixParameter(2,f_bw->GetParameter(5));

    //      h_mMass[KEY]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
    //    }
    //  }
    //}

    ////////////////////////////////////////////////////////////////////
    // Yields
    TH1FMap h_mYield_SE, h_mYield_ME, h_mYield;
    for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
    {
      string KEY_Yield_SE = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SE",i_cent,vmsa::mPID[pid].c_str());
      string Hist_Yield_SE = KEY_Yield_SE;
      h_mYield_SE[KEY_Yield_SE] = (TH1F*)File_SE->Get(Hist_Yield_SE.c_str())->Clone();
      Int_t Norm_bin_start = h_mYield_SE[KEY_Yield_SE]->FindBin(vmsa::Norm_Start[pid][0]);
      Int_t Norm_bin_stop  = h_mYield_SE[KEY_Yield_SE]->FindBin(vmsa::Norm_Stop[pid][0]);
      Float_t Inte_SE = h_mYield_SE[KEY_Yield_SE]->Integral(Norm_bin_start,Norm_bin_stop);

      string KEY_Yield_ME = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_ME",i_cent,vmsa::mPID[pid].c_str());
      string Hist_Yield_ME = KEY_Yield_ME;
      h_mYield_ME[KEY_Yield_ME] = (TH1F*)File_ME->Get(Hist_Yield_ME.c_str())->Clone();
      Float_t Inte_ME = h_mYield_ME[KEY_Yield_ME]->Integral(Norm_bin_start,Norm_bin_stop);
      h_mYield_ME[KEY_Yield_ME]->Scale(Inte_SE/Inte_ME);

      string KEY_Yield = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SM",i_cent,vmsa::mPID[pid].c_str());
      h_mYield[KEY_Yield] = (TH1F*)h_mYield_SE[KEY_Yield_SE]->Clone();
      h_mYield[KEY_Yield]->Add(h_mYield_ME[KEY_Yield_ME],-1.0); 
      cout << KEY_Yield_SE << endl;
    }


    TH1FMap h_mYield_SM; // for QA plot only
    vecFMap ParYield_SM;

    for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
    {
      string KEY_Yield = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SM",i_cent,vmsa::mPID[pid].c_str());
      h_mYield_SM[KEY_Yield] = (TH1F*)h_mYield[KEY_Yield]->Clone();
      TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      for(Int_t i_par = 0; i_par < 7; i_par++)
      {
        f_bw->ReleaseParameter(i_par);
      }
      f_bw->SetParameter(0,vmsa::InvMass[pid]);
      f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.01,vmsa::InvMass[pid]+0.01);
      f_bw->SetParameter(1,vmsa::Width[pid]);
      f_bw->SetParameter(2,10000);
      f_bw->SetParameter(3,6000);
      //f_bw->SetParameter(4,0.5);
      //f_bw->SetParameter(5,0.1);
      f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
      ParYield_SM[KEY_Yield].clear();
      h_mYield[KEY_Yield]->Fit(f_bw,"QNR");
      for(Int_t n_par = 0; n_par < 7; n_par++)
      {
        ParYield_SM[KEY_Yield].push_back(static_cast<Float_t>(f_bw->GetParameter(n_par)));
      }
    }

    vecFMap yields_Gaus, yields_BW;
    for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
    {
      string KEY_Yield = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SM",i_cent,vmsa::mPID[pid].c_str());
      TF1 *f_yields_bw = new TF1("f_yields_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      f_yields_bw->SetParameter(0,ParYield_SM[KEY_Yield][0]);
      f_yields_bw->SetParameter(1,ParYield_SM[KEY_Yield][1]);
      f_yields_bw->SetParameter(2,ParYield_SM[KEY_Yield][2]);
      f_yields_bw->SetParameter(3,ParYield_SM[KEY_Yield][3]);
      f_yields_bw->SetParameter(4,ParYield_SM[KEY_Yield][4]);
      f_yields_bw->SetParameter(5,ParYield_SM[KEY_Yield][5]);
      f_yields_bw->SetParameter(6,ParYield_SM[KEY_Yield][6]);
      f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
      TFitResultPtr result = h_mYield[KEY_Yield]->Fit(f_yields_bw,"MQNRIS");

      TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
      f_bg->SetParameter(0,ParYield_SM[KEY_Yield][3]);
      f_bg->SetParameter(1,ParYield_SM[KEY_Yield][4]);
      f_bg->SetParameter(2,ParYield_SM[KEY_Yield][5]);
      f_bg->SetParameter(3,ParYield_SM[KEY_Yield][6]);
      f_bg->SetParError(0,ParYield_SM[KEY_Yield][3]);
      f_bg->SetParError(1,ParYield_SM[KEY_Yield][4]);
      f_bg->SetParError(2,ParYield_SM[KEY_Yield][5]);
      f_bg->SetParError(3,ParYield_SM[KEY_Yield][6]);
      
      double params[4] = {result->GetParams()[3],result->GetParams()[4],result->GetParams()[5],result->GetParams()[6]};
      TMatrixDSym covArr(4);
      covArr(0,0) = result->GetCovarianceMatrix()(3,3);
      covArr(0,1) = result->GetCovarianceMatrix()(3,4);
      covArr(0,2) = result->GetCovarianceMatrix()(3,5);
      covArr(0,3) = result->GetCovarianceMatrix()(3,6);
      covArr(1,0) = result->GetCovarianceMatrix()(4,3);
      covArr(1,1) = result->GetCovarianceMatrix()(4,4);
      covArr(1,2) = result->GetCovarianceMatrix()(4,5);
      covArr(1,3) = result->GetCovarianceMatrix()(4,6);
      covArr(2,0) = result->GetCovarianceMatrix()(5,3);
      covArr(2,1) = result->GetCovarianceMatrix()(5,4);
      covArr(2,2) = result->GetCovarianceMatrix()(5,5);
      covArr(2,3) = result->GetCovarianceMatrix()(5,6);
      covArr(3,0) = result->GetCovarianceMatrix()(6,3);
      covArr(3,1) = result->GetCovarianceMatrix()(6,4);
      covArr(3,2) = result->GetCovarianceMatrix()(6,5);
      covArr(3,3) = result->GetCovarianceMatrix()(6,6);


      // counting for guassian
      Float_t counts_gaus = 0.0;
      Float_t errors_gaus = 0.0;
      Int_t bin_start = h_mYield[KEY_Yield]->FindBin(ParYield_SM[KEY_Yield][0]-vmsa::nSigVec*ParYield_SM[KEY_Yield][1]);
      Int_t bin_stop  = h_mYield[KEY_Yield]->FindBin(ParYield_SM[KEY_Yield][0]+vmsa::nSigVec*ParYield_SM[KEY_Yield][1]);
      for(Int_t i_bin = bin_start; i_bin <= bin_stop; i_bin++)
      {
        counts_gaus += h_mYield[KEY_Yield]->GetBinContent(i_bin);
        errors_gaus += h_mYield[KEY_Yield]->GetBinError(i_bin)*h_mYield[KEY_Yield]->GetBinError(i_bin);
      }

      // integrating for breit wigner
      Float_t bin_width = h_mYield[KEY_Yield]->GetBinWidth(1);
      Float_t Inte_start = ParYield_SM[KEY_Yield][0]-vmsa::nSigVec*ParYield_SM[KEY_Yield][1]-0.5*bin_width;
      Float_t Inte_stop  = ParYield_SM[KEY_Yield][0]+vmsa::nSigVec*ParYield_SM[KEY_Yield][1]+0.5*bin_width;
      Float_t counts_bw = f_yields_bw->Integral(Inte_start,Inte_stop)/bin_width;
      Float_t errors_bw = f_yields_bw->IntegralError(Inte_start,Inte_stop)/bin_width;

      float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
      float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width;

      yields_Gaus[KEY_Yield].clear();
      yields_Gaus[KEY_Yield].push_back(static_cast<Float_t>(counts_gaus-counts_bg));
      yields_Gaus[KEY_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors_gaus+errors_bg*errors_bg)));
      yields_BW[KEY_Yield].clear();
      yields_BW[KEY_Yield].push_back(static_cast<Float_t>(counts_bw-counts_bg));
      yields_BW[KEY_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg)));
    }

    // QA: different counting method: bin counting vs breit wigner integrating
    TCanvas *c_Yields_counts = new TCanvas("c_Yields_counts","c_Yields_counts",10,10,900,900);
    c_Yields_counts->Divide(3,3);
    for(Int_t i_cent = 0; i_cent < 9; i_cent++)
    {
      c_Yields_counts->cd(i_cent+1);
      c_Yields_counts->cd(i_cent+1)->SetLeftMargin(0.20);
      c_Yields_counts->cd(i_cent+1)->SetBottomMargin(0.20);
      c_Yields_counts->cd(i_cent+1)->SetTicks(1,1);
      c_Yields_counts->cd(i_cent+1)->SetGrid(0,0);
  
      string KEY_Yield = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SM",i_cent,vmsa::mPID[pid].c_str());
      h_mYield[KEY_Yield]->SetTitle("");
      h_mYield[KEY_Yield]->SetStats(0);
      h_mYield[KEY_Yield]->SetMarkerStyle(24);
      h_mYield[KEY_Yield]->SetMarkerColor(kGray+3);
      h_mYield[KEY_Yield]->SetMarkerSize(0.8);
      h_mYield[KEY_Yield]->GetXaxis()->SetNdivisions(505,'N');
      h_mYield[KEY_Yield]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
      h_mYield[KEY_Yield]->GetXaxis()->CenterTitle();
      h_mYield[KEY_Yield]->GetYaxis()->SetTitle("Counts");
      h_mYield[KEY_Yield]->GetYaxis()->CenterTitle();
      h_mYield[KEY_Yield]->GetYaxis()->SetTitleOffset(1.2);
      h_mYield[KEY_Yield]->DrawCopy("pE");
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  
      TF1 *f_yields_bw = new TF1("f_yields_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      f_yields_bw->SetParameter(0,ParYield_SM[KEY_Yield][0]);
      f_yields_bw->SetParameter(1,ParYield_SM[KEY_Yield][1]);
      f_yields_bw->SetParameter(2,ParYield_SM[KEY_Yield][2]);
      f_yields_bw->SetParameter(3,ParYield_SM[KEY_Yield][3]);
      f_yields_bw->SetParameter(4,ParYield_SM[KEY_Yield][4]);
      f_yields_bw->SetParameter(5,ParYield_SM[KEY_Yield][5]);
      f_yields_bw->SetParameter(6,ParYield_SM[KEY_Yield][6]);
      f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
      f_yields_bw->SetLineColor(2);
      f_yields_bw->SetLineStyle(1);
      f_yields_bw->SetLineWidth(2);
      f_yields_bw->DrawCopy("l same");
  
      Float_t x1 = ParYield_SM[KEY_Yield][0] - vmsa::nSigVec*ParYield_SM[KEY_Yield][1];
      Float_t x2 = ParYield_SM[KEY_Yield][0] + vmsa::nSigVec*ParYield_SM[KEY_Yield][1];
      Float_t y = h_mYield[KEY_Yield]->GetBinContent(h_mYield[KEY_Yield]->FindBin(ParYield_SM[KEY_Yield][0]));
      PlotLine(x1,x1,0,y,4,2,2);
      PlotLine(x2,x2,0,y,4,2,2);
    }

    Float_t yields_total_gaus[4] = {0.0};
    Float_t yields_total_bw[4]   = {0.0};
    // calculate final resolution correction factors and correct flow
    for(Int_t i_cent = vmsa::Cent_start; i_cent < 9; i_cent++) // centrality bin loop
    { 
      string KEY_Yield = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SM",i_cent,vmsa::mPID[pid].c_str());

      if(i_cent >= 0 && i_cent <= 8) // calculate resolution and total yields in selected centrality bin
      {
        yields_total_gaus[0] += yields_Gaus[KEY_Yield][0];
        yields_total_bw[0] += yields_BW[KEY_Yield][0];
      }
      if(i_cent == 8 || i_cent == 7) // calculate resolution and total yields in selected centrality bin
      {
        yields_total_gaus[1] += yields_Gaus[KEY_Yield][0];
        yields_total_bw[1] += yields_BW[KEY_Yield][0];
      }
      if(i_cent <= 6 && i_cent >= 4) // calculate resolution and total yields in selected centrality bin
      {
        yields_total_gaus[2] += yields_Gaus[KEY_Yield][0];
        yields_total_bw[2] += yields_BW[KEY_Yield][0];
      }
      if(i_cent <= 3 && i_cent >= 0) // calculate resolution and total yields in selected centrality bin
      {
        yields_total_gaus[3] += yields_Gaus[KEY_Yield][0];
        yields_total_bw[3] += yields_BW[KEY_Yield][0];
      }
    }

    Float_t mean_res_gaus[4] = {0.0};
    Float_t mean_res_bw[4] = {0.0};
    for(Int_t i_cent = vmsa::Cent_start; i_cent < 9; i_cent++) // centrality bin loop
    {   

      string KEY_Yield = Form("Yields_Centrality_%d_Dca_0_Sig_0_NHit_0_%s_SM",i_cent,vmsa::mPID[pid].c_str());

      if(i_cent >= 0 && i_cent <= 8) // calculate resolution and total yields in selected centrality bin
      {
        mean_res_gaus[0] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[0]);
        mean_res_bw[0] += yields_BW[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_bw[0]);
      }
      if(i_cent == 8 || i_cent == 7) // calculate resolution and total yields in selected centrality bin
      {
        mean_res_gaus[1] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[1]);
        mean_res_bw[1] += yields_BW[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_bw[1]);
      }
      if(i_cent <= 6 && i_cent >= 4) // calculate resolution and total yields in selected centrality bin
      {
        mean_res_gaus[2] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[2]);
        mean_res_bw[2] += yields_BW[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_bw[2]);
      }
      if(i_cent <= 3 && i_cent >= 0) // calculate resolution and total yields in selected centrality bin
      {
        mean_res_gaus[3] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[3]);
        mean_res_bw[3] += yields_BW[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_bw[3]);
      }
    }

    for(Int_t i_cent = vmsa::Cent_start; i_cent < 4; i_cent++) // centrality bin loop
    {   
      cout << "centrality_bin = " << i_cent << ", mean_res_gaus = " << mean_res_gaus[i_cent] << endl;
      cout << "centrality_bin = " << i_cent << ", mean_res_bw = " << mean_res_bw[i_cent] << endl;
    }
      


    ///////////////////////////////////////////////////////////

#if _PlotQA_
    // QA plots for phi differential InvMass after linear background subtraction
    TCanvas *c_mMass_phi_diff = new TCanvas("c_mMass_phi_diff","c_mMass_phi_diff",10,10,2500,2500);
    c_mMass_phi_diff->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,9,0,vmsa::mPID[pid].c_str());
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("PE");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }
    c_mMass_phi_diff->SaveAs("./figures/c_mMass_phi_diff.pdf");
#endif
  //}

  TH1FMap h_mMass_total; // cos(theta*) integrated InvMass after linear background subtraction for bw fits to extract yields 
  vecFMap ParBW;

  for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
    {
      string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
      for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
      {
        string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
        if(i_theta == 0) h_mMass_total[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone();
        else h_mMass_total[KEY_theta]->Add(h_mMass[KEY],1.0);
      }
      TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      f_bw->SetParameter(0,vmsa::InvMass[pid]);
      f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.01,vmsa::InvMass[pid]+0.01);
      f_bw->SetParameter(1,vmsa::Width[pid]);
      f_bw->SetParameter(2,10000);
      f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
      h_mMass_total[KEY_theta]->Fit(f_bw,"MQNR");
      ParBW[KEY_theta].clear();
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(0)));
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(1)));
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(2)));
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(3)));
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(4)));
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(5)));
      ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(6)));
    }
  }

#if _PlotQA_
  // QA: bw fits to phi integrated InvMass
  TCanvas *c_mMass_bw = new TCanvas("c_mMass_bw","c_mMass_bw",10,10,1400,1400);
  c_mMass_bw->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
  {
    c_mMass_bw->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
    c_mMass_bw->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_bw->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_bw->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
    c_mMass_bw->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,9,vmsa::mPID[pid].c_str());
    h_mMass_total[KEY_theta_QA]->SetMarkerColor(1);
    h_mMass_total[KEY_theta_QA]->SetMarkerStyle(24);
    h_mMass_total[KEY_theta_QA]->SetMarkerSize(0.8);
    h_mMass_total[KEY_theta_QA]->DrawCopy("PE");
    TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
    f_bw->SetParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->SetParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]);
    f_bw->SetParameter(3,ParBW[KEY_theta_QA][3]);
    f_bw->SetParameter(4,ParBW[KEY_theta_QA][4]);
    f_bw->SetParameter(5,ParBW[KEY_theta_QA][5]);
    f_bw->SetParameter(6,ParBW[KEY_theta_QA][6]);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }
#endif

  double resMean = 0.0;
  double resWeightI[2][7][10];
  double resWeight[2][7];
  double phiMesons[2][7][10];
  memset(phiMesons, 0.0, sizeof phiMesons);

  // calculate counts and errors for cos(theta*) bin with bin counting and integrating
  TH1FMap h_mCounts, h_mv2, h_mv2_corr;
  vecFMap ParSpin_Count, ParSpin_BW;
  for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
  {
    string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
    h_mv2[KEY_v2_Count] = new TH1F(KEY_v2_Count.c_str(),KEY_v2_Count.c_str(),100,-0.05,9.95);
    string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
    h_mv2[KEY_v2_BW] = new TH1F(KEY_v2_BW.c_str(),KEY_v2_BW.c_str(),100,-0.05,9.95);

    string KEY_v2_Count_corr = Form("v2_Corr_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
    h_mv2_corr[KEY_v2_Count_corr] = new TH1F(KEY_v2_Count_corr.c_str(),KEY_v2_Count_corr.c_str(),100,-0.05,9.95);
    string KEY_v2_BW_corr = Form("v2_Corr_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
    h_mv2_corr[KEY_v2_BW_corr] = new TH1F(KEY_v2_BW_corr.c_str(),KEY_v2_BW_corr.c_str(),100,-0.05,9.95);

    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
      h_mCounts[KEY_Count] = new TH1F(KEY_Count.c_str(),KEY_Count.c_str(),10,0.0,TMath::Pi()/2.0);
      string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
      h_mCounts[KEY_BW] = new TH1F(KEY_BW.c_str(),KEY_BW.c_str(),10,0.0,TMath::Pi()/2.0);
      string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
      for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
      {
        // bin counting
        string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
        float counts = 0.0;
        float errors = 0.0;
        float bin_center = (vmsa::PhiPsi_low[i_theta] + vmsa::PhiPsi_up[i_theta])/2.0;
        int bin_start = h_mMass[KEY]->FindBin(ParBW[KEY_theta][0]-vmsa::nSigVec*ParBW[KEY_theta][1]);
        int bin_stop  = h_mMass[KEY]->FindBin(ParBW[KEY_theta][0]+vmsa::nSigVec*ParBW[KEY_theta][1]);
        for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
        {
          counts += h_mMass[KEY]->GetBinContent(i_bin);
          errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
        }
      
        // breit wigner integrating
        TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
        f_bw->FixParameter(0,ParBW[KEY_theta][0]);
        f_bw->FixParameter(1,ParBW[KEY_theta][1]);
        f_bw->SetParameter(2,ParBW[KEY_theta][2]/10.0);
        f_bw->FixParameter(3,ParBW[KEY_theta][3]/10.0);
        f_bw->FixParameter(4,ParBW[KEY_theta][4]/10.0);
        f_bw->SetParameter(5,ParBW[KEY_theta][5]/10.0);
        f_bw->SetParameter(6,ParBW[KEY_theta][6]/10.0);
        f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
        TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"NMQRIS");

        TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);;
        f_bg->SetParameter(0,f_bw->GetParameter(3));
        f_bg->SetParameter(1,f_bw->GetParameter(4));
        f_bg->SetParameter(2,f_bw->GetParameter(5));
        f_bg->SetParameter(3,f_bw->GetParameter(6));
        f_bg->SetParError(0,f_bw->GetParError(3));
        f_bg->SetParError(1,f_bw->GetParError(4));
        f_bg->SetParError(2,f_bw->GetParError(5));
        f_bg->SetParError(3,f_bw->GetParError(6));

        double params[4] = {result->GetParams()[3],result->GetParams()[4],result->GetParams()[5],result->GetParams()[6]};
        TMatrixDSym covArr(4);
        covArr(0,0) = result->GetCovarianceMatrix()(3,3);
        covArr(0,1) = result->GetCovarianceMatrix()(3,4);
        covArr(0,2) = result->GetCovarianceMatrix()(3,5);
        covArr(0,3) = result->GetCovarianceMatrix()(3,6);
        covArr(1,0) = result->GetCovarianceMatrix()(4,3);
        covArr(1,1) = result->GetCovarianceMatrix()(4,4);
        covArr(1,2) = result->GetCovarianceMatrix()(4,5);
        covArr(1,3) = result->GetCovarianceMatrix()(4,6);
        covArr(2,0) = result->GetCovarianceMatrix()(5,3);
        covArr(2,1) = result->GetCovarianceMatrix()(5,4);
        covArr(2,2) = result->GetCovarianceMatrix()(5,5);
        covArr(2,3) = result->GetCovarianceMatrix()(5,6);
        covArr(3,0) = result->GetCovarianceMatrix()(6,3);
        covArr(3,1) = result->GetCovarianceMatrix()(6,4);
        covArr(3,2) = result->GetCovarianceMatrix()(6,5);
        covArr(3,3) = result->GetCovarianceMatrix()(6,6);


        float bin_width = h_mMass[KEY]->GetBinWidth(1);
        float Inte_start = ParBW[KEY_theta][0]-vmsa::nSigVec*ParBW[KEY_theta][1]-0.5*bin_width;
        float Inte_stop  = ParBW[KEY_theta][0]+vmsa::nSigVec*ParBW[KEY_theta][1]+0.5*bin_width;
        float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
        float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
     
        float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
        float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width;
        h_mCounts[KEY_Count]->SetBinContent(h_mCounts[KEY_Count]->FindBin(bin_center),counts-counts_bg);
        h_mCounts[KEY_Count]->SetBinError(h_mCounts[KEY_Count]->FindBin(bin_center),TMath::Sqrt(errors+errors_bg*errors_bg));
        h_mCounts[KEY_BW]->SetBinContent(h_mCounts[KEY_BW]->FindBin(bin_center),counts_bw-counts_bg);
        h_mCounts[KEY_BW]->SetBinError(h_mCounts[KEY_BW]->FindBin(bin_center),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));  
      }
      float pt_mean = (vmsa::pt_lowKSv2[energy][i_pt]+vmsa::pt_upKSv2[energy][i_pt])/2.0;
      
      TF1 *f_cos_count = new TF1("f_cos_count",flow,0.0,TMath::Pi()/2.0,2);
      f_cos_count->SetParameter(0,0.33);
      f_cos_count->SetParameter(1,1.0);
      h_mCounts[KEY_Count]->Fit(f_cos_count,"NQMRI");
      ParSpin_Count[KEY_Count].clear();
      ParSpin_Count[KEY_Count].push_back(static_cast<float>(f_cos_count->GetParameter(0)));
      ParSpin_Count[KEY_Count].push_back(static_cast<float>(f_cos_count->GetParameter(1)));
      h_mv2[KEY_v2_Count]->SetBinContent(h_mv2[KEY_v2_Count]->FindBin(pt_mean),f_cos_count->GetParameter(0));
      h_mv2[KEY_v2_Count]->SetBinError(h_mv2[KEY_v2_Count]->FindBin(pt_mean),f_cos_count->GetParError(0));
      
      TF1 *f_cos_bw = new TF1("f_cos_bw",flow,0.0,TMath::Pi()/2.0,2);
      f_cos_bw->SetParameter(0,0.33);
      f_cos_bw->SetParameter(1,1.0);
      h_mCounts[KEY_BW]->Fit(f_cos_bw,"NQMI");
      ParSpin_BW[KEY_BW].clear();
      ParSpin_BW[KEY_BW].push_back(static_cast<float>(f_cos_bw->GetParameter(0)));
      ParSpin_BW[KEY_BW].push_back(static_cast<float>(f_cos_bw->GetParameter(1)));
      h_mv2[KEY_v2_BW]->SetBinContent(h_mv2[KEY_v2_BW]->FindBin(pt_mean),f_cos_bw->GetParameter(0));
      h_mv2[KEY_v2_BW]->SetBinError(h_mv2[KEY_v2_BW]->FindBin(pt_mean),f_cos_bw->GetParError(0));
      
      if(i_cent > 8)
      {
        cout << "cent     = " << i_cent << endl;
        cout << "count v2 = " << f_cos_count->GetParameter(0)*mean_res_gaus[i_cent-9] << endl;
        cout << "   bw v2 = " << f_cos_bw->GetParameter(0)*mean_res_bw[i_cent-9] << endl;
        h_mv2_corr[KEY_v2_Count_corr]->SetBinContent(h_mv2_corr[KEY_v2_Count_corr]->FindBin(pt_mean),f_cos_count->GetParameter(0)*mean_res_gaus[i_cent-9]);
        h_mv2_corr[KEY_v2_Count_corr]->SetBinError(h_mv2_corr[KEY_v2_Count_corr]->FindBin(pt_mean),f_cos_count->GetParError(0)*mean_res_gaus[i_cent-9]);
        h_mv2_corr[KEY_v2_BW_corr]->SetBinContent(h_mv2_corr[KEY_v2_BW_corr]->FindBin(pt_mean),f_cos_bw->GetParameter(0)*mean_res_bw[i_cent-9]);
        h_mv2_corr[KEY_v2_BW_corr]->SetBinError(h_mv2_corr[KEY_v2_BW_corr]->FindBin(pt_mean),f_cos_bw->GetParError(0)*mean_res_bw[i_cent-9]);
      }
    }
  }

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",2,9,vmsa::mPID[pid].c_str());
  TCanvas *c_mMass_psi = new TCanvas("c_mMass_psi","c_mMass_psi",10,10,800,800);
  string Title_PhiPsi[3] = {"2/7 < cos(#theta*) < 3/7","3/7 < cos(#theta*) < 4/7","4/7 < cos(#theta*) < 5/7"};
  c_mMass_psi->Divide(2,2);
  for(int i_theta = 2; i_theta < 5; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi->cd(i_theta-1)->SetLeftMargin(0.15);
    c_mMass_psi->cd(i_theta-1)->SetBottomMargin(0.15);
    c_mMass_psi->cd(i_theta-1)->SetTicks(1,1);
    c_mMass_psi->cd(i_theta-1)->SetGrid(0,0);
    c_mMass_psi->cd(i_theta-1);
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",2,9,i_theta,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetTitle("");
    h_mMass[KEY_QA]->SetStats(0);
    h_mMass[KEY_QA]->GetXaxis()->SetRangeUser(0.0,5.0);
    h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
    h_mMass[KEY_QA]->GetXaxis()->SetTitle("InvMass(K^{+},K^{-}) (GeV/c^{2})");
    h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
    h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.05);
    h_mMass[KEY_QA]->GetXaxis()->CenterTitle();

    h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(0.0,1.1*h_mMass[KEY_QA]->GetMaximum());
    h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
    h_mMass[KEY_QA]->GetYaxis()->SetTitle("Counts");
    h_mMass[KEY_QA]->GetYaxis()->SetTitleOffset(1.2);
    h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.04);
    h_mMass[KEY_QA]->GetYaxis()->CenterTitle();
    h_mMass[KEY_QA]->SetMarkerColor(1);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(0.8);
    h_mMass[KEY_QA]->DrawCopy("hE");

    TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
    f_bw->FixParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->FixParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]/10.0);
    f_bw->SetParameter(3,ParBW[KEY_theta_QA][3]/10.0);
    f_bw->SetParameter(4,ParBW[KEY_theta_QA][4]/10.0);
    f_bw->SetParameter(5,ParBW[KEY_theta_QA][5]/10.0);
    f_bw->SetParameter(6,ParBW[KEY_theta_QA][6]/10.0);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    h_mMass[KEY_QA]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    float x1 = ParBW[KEY_theta_QA][0] - vmsa::nSigVec*ParBW[KEY_theta_QA][1];
    float x2 = ParBW[KEY_theta_QA][0] + vmsa::nSigVec*ParBW[KEY_theta_QA][1];
    float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(ParBW[KEY_theta_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }

  c_mMass_psi->cd(4);
  c_mMass_psi->cd(4)->SetLeftMargin(0.15);
  c_mMass_psi->cd(4)->SetBottomMargin(0.15);
  c_mMass_psi->cd(4)->SetTicks(1,1);
  c_mMass_psi->cd(4)->SetGrid(0,0);
  string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",2,9,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_Count]->SetStats(0);
  // h_mCounts[KEY_Count]->SetTitle("20-60%");
  h_mCounts[KEY_Count]->SetTitle("");
  h_mCounts[KEY_Count]->SetTitleSize(0.08);
  h_mCounts[KEY_Count]->GetXaxis()->SetTitle("cos(#theta*) (w.r.t. 2^{nd} event plane)");
  h_mCounts[KEY_Count]->GetXaxis()->SetTitleSize(0.06);
  h_mCounts[KEY_Count]->GetXaxis()->CenterTitle();
  h_mCounts[KEY_Count]->GetXaxis()->SetNdivisions(505);
  h_mCounts[KEY_Count]->GetXaxis()->SetLabelSize(0.05);

  h_mCounts[KEY_Count]->GetYaxis()->SetTitle("Yields");
  h_mCounts[KEY_Count]->GetYaxis()->SetTitleSize(0.06);
  h_mCounts[KEY_Count]->GetYaxis()->CenterTitle();
  h_mCounts[KEY_Count]->GetYaxis()->SetNdivisions(505);
  h_mCounts[KEY_Count]->GetYaxis()->SetLabelSize(0.05);
  h_mCounts[KEY_Count]->SetLineColor(4);
  h_mCounts[KEY_Count]->SetMarkerColor(4);
  h_mCounts[KEY_Count]->SetMarkerStyle(24);
  h_mCounts[KEY_Count]->SetMarkerSize(1.2);
  h_mCounts[KEY_Count]->DrawCopy("pE");
  TF1 *f_cos_count_QA = new TF1("f_cos_count_QA",flow,0.0,TMath::Pi()/2.0,2);
  f_cos_count_QA->FixParameter(0,ParSpin_Count[KEY_Count][0]);
  f_cos_count_QA->FixParameter(1,ParSpin_Count[KEY_Count][1]);
  f_cos_count_QA->SetLineStyle(2);
  f_cos_count_QA->SetLineColor(4);
  f_cos_count_QA->DrawCopy("l same");

  string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",2,9,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(1.2);
  h_mCounts[KEY_BW]->DrawCopy("pE same");
  TF1 *f_cos_bw_QA = new TF1("f_cos_bw_QA",flow,0.0,TMath::Pi()/2.0,2);
  f_cos_bw_QA->FixParameter(0,ParSpin_BW[KEY_BW][0]);
  f_cos_bw_QA->FixParameter(1,ParSpin_BW[KEY_BW][1]);
  f_cos_bw_QA->SetLineStyle(2);
  f_cos_bw_QA->SetLineColor(2);
  f_cos_bw_QA->DrawCopy("l same");

  TLegend *leg_temp = new TLegend(0.18,0.3,0.55,0.5);
  leg_temp->SetFillColor(10);
  leg_temp->SetBorderSize(0.0);
  leg_temp->AddEntry(h_mCounts[KEY_Count],"bin counting","p");
  leg_temp->AddEntry(h_mCounts[KEY_BW],"Breit-Wigner","p");
  leg_temp->Draw("same");
  c_mMass_psi->SaveAs("./figures/phi_SpinAlighment.png");
#endif

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  TCanvas *c_mMass_psi_2 = new TCanvas("c_mMass_psi_2","c_mMass_psi_2",10,10,1000,1400);
  c_mMass_psi_2->Divide(3,4);
  for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi_2->cd(i_theta+1)->SetLeftMargin(0.15);
    c_mMass_psi_2->cd(i_theta+1)->SetBottomMargin(0.15);
    c_mMass_psi_2->cd(i_theta+1)->SetTicks(1,1);
    c_mMass_psi_2->cd(i_theta+1)->SetGrid(0,0);
    c_mMass_psi_2->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",2,9,i_theta,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetTitle(Form("#phi-Psi_{2} bin %d",i_theta));
    h_mMass[KEY_QA]->SetMarkerColor(1);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(0.8);
    h_mMass[KEY_QA]->DrawCopy("pE");

    TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
    f_bw->FixParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->FixParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]/10.0);
    f_bw->SetParameter(3,ParBW[KEY_theta_QA][3]/10.0);
    f_bw->SetParameter(4,ParBW[KEY_theta_QA][4]/10.0);
    f_bw->SetParameter(5,ParBW[KEY_theta_QA][5]/10.0);
    f_bw->SetParameter(6,ParBW[KEY_theta_QA][6]/10.0);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    h_mMass[KEY_QA]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    float x1 = ParBW[KEY_theta_QA][0] - vmsa::nSigVec*ParBW[KEY_theta_QA][1];
    float x2 = ParBW[KEY_theta_QA][0] + vmsa::nSigVec*ParBW[KEY_theta_QA][1];
    float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(ParBW[KEY_theta_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }

  c_mMass_psi_2->cd(11);
  h_mMass_total[KEY_theta_QA]->SetTitle("Integrated Yields");
  h_mMass_total[KEY_theta_QA]->SetMarkerColor(1);
  h_mMass_total[KEY_theta_QA]->SetMarkerStyle(24);
  h_mMass_total[KEY_theta_QA]->SetMarkerSize(0.8);
  h_mMass_total[KEY_theta_QA]->DrawCopy("pE");
  TF1 *f_bw_QA = new TF1("f_bw_QA",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
  f_bw_QA->FixParameter(0,ParBW[KEY_theta_QA][0]);
  f_bw_QA->FixParameter(1,ParBW[KEY_theta_QA][1]);
  f_bw_QA->FixParameter(2,ParBW[KEY_theta_QA][2]);
  f_bw_QA->FixParameter(3,ParBW[KEY_theta_QA][3]);
  f_bw_QA->FixParameter(4,ParBW[KEY_theta_QA][4]);
  f_bw_QA->FixParameter(5,ParBW[KEY_theta_QA][5]);
  f_bw_QA->FixParameter(6,ParBW[KEY_theta_QA][6]);
  f_bw_QA->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
  f_bw_QA->SetLineColor(2);
  f_bw_QA->Draw("l same");
  string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][2],vmsa::pt_upKSv2[energy][2]);
  plotTopLegend((char*)pT_range.c_str(),0.15,0.7,0.08,1,0.0,42,1);
  PlotLine(0.98,1.05,0.0,0.0,1,2,2);


  c_mMass_psi_2->cd(12);
  h_mCounts[KEY_Count]->SetStats(0);
  h_mCounts[KEY_Count]->SetTitle("20-60%");
  h_mCounts[KEY_Count]->SetTitleSize(0.08);
  h_mCounts[KEY_Count]->SetLineColor(4);
  h_mCounts[KEY_Count]->SetMarkerColor(4);
  h_mCounts[KEY_Count]->SetMarkerStyle(24);
  h_mCounts[KEY_Count]->SetMarkerSize(0.8);
  h_mCounts[KEY_Count]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  h_mCounts[KEY_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[KEY_Count]->GetXaxis()->CenterTitle();
  h_mCounts[KEY_Count]->DrawCopy("pE");
  TF1 *f_cos_count_QA_2 = new TF1("f_cos_count_QA_2",flow,0.0,TMath::Pi()/2.0,2);
  f_cos_count_QA_2->FixParameter(0,ParSpin_Count[KEY_Count][0]);
  f_cos_count_QA_2->FixParameter(1,ParSpin_Count[KEY_Count][1]);
  f_cos_count_QA_2->SetLineStyle(2);
  f_cos_count_QA_2->SetLineColor(4);
  f_cos_count_QA_2->DrawCopy("l same");

  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(0.8);
  h_mCounts[KEY_BW]->DrawCopy("pE same");
  TF1 *f_cos_bw_QA_2 = new TF1("f_cos_bw_QA_2",flow,0.0,TMath::Pi()/2.0,2);
  f_cos_bw_QA_2->FixParameter(0,ParSpin_BW[KEY_BW][0]);
  f_cos_bw_QA_2->FixParameter(1,ParSpin_BW[KEY_BW][1]);
  f_cos_bw_QA_2->SetLineStyle(2);
  f_cos_bw_QA_2->SetLineColor(2);
  f_cos_bw_QA_2->DrawCopy("l same");

  TLegend *leg_temp_2 = new TLegend(0.5,0.6,0.8,0.8);
  leg_temp_2->SetFillColor(10);
  leg_temp_2->SetBorderSize(0.0);
  leg_temp_2->AddEntry(h_mCounts[KEY_Count],"bin counting","p");
  leg_temp_2->AddEntry(h_mCounts[KEY_BW],"breit wigner","p");
  leg_temp_2->Draw("same");
  c_mMass_psi_2->SaveAs("./figures/phi_SpinAlighment.pdf");
#endif

  // set final pt and rho00 to one TGraphAsymmErrors
  double weightedMeanC = 0.0;
  double weightedErrorC = 0.0;
  double weightedMeanI = 0.0;
  double weightedErrorI = 0.0;


  TGraMap g_mv2;
  TGraMap g_mv2_corr;
  for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
  {
    string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
    g_mv2[KEY_v2_Count] = new TGraphAsymmErrors();
    string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
    g_mv2[KEY_v2_BW] = new TGraphAsymmErrors();

    string KEY_v2_Count_corr;  
    string KEY_v2_BW_corr;
    if(i_cent > 8)
    {          
      KEY_v2_Count_corr = Form("v2_Corr_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      g_mv2_corr[KEY_v2_Count_corr] = new TGraphAsymmErrors();
      KEY_v2_BW_corr = Form("v2_Corr_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
      g_mv2_corr[KEY_v2_BW_corr] = new TGraphAsymmErrors();
    }

    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++)
    {
 
      if(i_cent <= 8)
      {

        float pt_mean = (vmsa::pt_lowKSv2[energy][i_pt]+vmsa::pt_upKSv2[energy][i_pt])/2.0;

        float count_content = h_mv2[KEY_v2_Count]->GetBinContent(h_mv2[KEY_v2_Count]->FindBin(pt_mean)); // bin counting
        float count_error   = h_mv2[KEY_v2_Count]->GetBinError(h_mv2[KEY_v2_Count]->FindBin(pt_mean));
 
        float bw_content = h_mv2[KEY_v2_BW]->GetBinContent(h_mv2[KEY_v2_BW]->FindBin(pt_mean)); // breit wigner fits
        float bw_error   = h_mv2[KEY_v2_BW]->GetBinError(h_mv2[KEY_v2_BW]->FindBin(pt_mean));

        g_mv2[KEY_v2_Count]->SetPoint(i_pt,pt_mean+0.05,count_content);
        g_mv2[KEY_v2_Count]->SetPointError(i_pt,0.0,0.0,count_error,count_error);
        g_mv2[KEY_v2_BW]->SetPoint(i_pt,pt_mean,bw_content);
        g_mv2[KEY_v2_BW]->SetPointError(i_pt,0.0,0.0,bw_error,bw_error);
 
        string Name_count = Form("v2_%s_Centrality_%d_Count",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
        g_mv2[KEY_v2_Count]->SetName(Name_count.c_str());
        g_mv2[KEY_v2_Count]->SetMarkerStyle(24);
        g_mv2[KEY_v2_Count]->SetMarkerColor(4);

        string Name_bw = Form("v2_%s_Centrality_%d_BW",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
        g_mv2[KEY_v2_BW]->SetName(Name_bw.c_str());
        g_mv2[KEY_v2_BW]->SetMarkerStyle(24);
        g_mv2[KEY_v2_BW]->SetMarkerColor(2);
      }
      if(i_cent > 8)
      {  
        float pt_mean = (vmsa::pt_lowKSv2[energy][i_pt]+vmsa::pt_upKSv2[energy][i_pt])/2.0;

        float count_content = h_mv2_corr[KEY_v2_Count_corr]->GetBinContent(h_mv2_corr[KEY_v2_Count_corr]->FindBin(pt_mean)); // bin counting
        float count_error   = h_mv2_corr[KEY_v2_Count_corr]->GetBinError(h_mv2_corr[KEY_v2_Count_corr]->FindBin(pt_mean));

        float bw_content = h_mv2_corr[KEY_v2_BW_corr]->GetBinContent(h_mv2_corr[KEY_v2_BW_corr]->FindBin(pt_mean)); // breit wigner fits
        float bw_error   = h_mv2_corr[KEY_v2_BW_corr]->GetBinError(h_mv2_corr[KEY_v2_BW_corr]->FindBin(pt_mean));

        cout << "cent = " << i_cent << endl;
        cout << "count_content = " << count_content << endl;
        cout << "bw_content = " << bw_content << endl;
    
        g_mv2_corr[KEY_v2_Count_corr]->SetPoint(i_pt,pt_mean+0.05,count_content);
        g_mv2_corr[KEY_v2_Count_corr]->SetPointError(i_pt,0.0,0.0,count_error,count_error);
        g_mv2_corr[KEY_v2_BW_corr]->SetPoint(i_pt,pt_mean,bw_content);
        g_mv2_corr[KEY_v2_BW_corr]->SetPointError(i_pt,0.0,0.0,bw_error,bw_error);

        string Name_count = Form("v2_corr_%s_Centrality_%d_Count",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
        g_mv2_corr[KEY_v2_Count_corr]->SetName(Name_count.c_str());
        g_mv2_corr[KEY_v2_Count_corr]->SetMarkerStyle(24);
        g_mv2_corr[KEY_v2_Count_corr]->SetMarkerColor(4);

        string Name_bw = Form("v2_corr_%s_Centrality_%d_BW",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
        g_mv2_corr[KEY_v2_BW_corr]->SetName(Name_bw.c_str());
        g_mv2_corr[KEY_v2_BW_corr]->SetMarkerStyle(24);
        g_mv2_corr[KEY_v2_BW_corr]->SetMarkerColor(2);
      }        
    }
  }


  //TFile *besi_cd = TFile::Open("../data/BESI/HEPData-ins1395151-v2-root.root");
  //TDirectory *dirAll = (TDirectory*) besi_cd->Get("Table 337;1");
  //dirAll->cd();
  //TGraphAsymmErrors *besi19_All = (TGraphAsymmErrors*)dirAll->Get("Graph1D_y1;1");
  //TDirectory *dir0 = (TDirectory*) besi_cd->Get("Table 43;1");
  //dir0->cd();
  //TGraphAsymmErrors *besi19_0 = (TGraphAsymmErrors*)dir0->Get("Graph1D_y1;1");
  //TDirectory *dir10 = (TDirectory*) besi_cd->Get("Table 141;1");
  //dir10->cd();
  //TGraphAsymmErrors *besi19_10 = (TGraphAsymmErrors*)dir10->Get("Graph1D_y1;1");
  //TDirectory *dir40 = (TDirectory*) besi_cd->Get("Table 239;1");
  //dir40->cd();
  //TGraphAsymmErrors *besi19_40 = (TGraphAsymmErrors*)dir40->Get("Graph1D_y1;1");

  // QA Plots for rho00 vs. pt of bin counting and breit wigner integrating
  {//  0-80%
  TCanvas *c_v2_pT = new TCanvas("c_v2_pT_0","c_v2_pT_0",10,10,800,800);
  c_v2_pT->cd();
  c_v2_pT->cd()->SetLeftMargin(0.15);
  c_v2_pT->cd()->SetBottomMargin(0.15);
  c_v2_pT->cd()->SetTicks(1,1);
  c_v2_pT->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(0.0,0.3);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");
  //PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  string KEY_v2_Count_QA = Form("v2_Corr_Count_Centrality_%d_%s",9,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_Count_QA] ,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.28,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.275,0.03,1,0.0,42,0);

  string KEY_v2_BW_QA = Form("v2_Corr_BW_Centrality_%d_%s",9,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_BW_QA] ,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.255,0.03,1,0.0,42,0);

  //string leg_besi = "BES-I";
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19_All,25,1,1.1);
  //Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,25,1,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.235,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.275,0.04,1,0.0,42,0);
  string leg_centrality = "0%-80%";
  plotTopLegend((char*)leg_centrality.c_str(),3.1,0.245,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/v2_pT_%s_0-80.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_v2_pT->SaveAs(figure_name.c_str());
  }
  {
  TCanvas *c_v2_pT = new TCanvas("c_v2_pT_1","c_v2_pT_1",10,10,800,800);
  c_v2_pT->cd();
  c_v2_pT->cd()->SetLeftMargin(0.15);
  c_v2_pT->cd()->SetBottomMargin(0.15);
  c_v2_pT->cd()->SetTicks(1,1);
  c_v2_pT->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(0.0,0.3);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");

  string KEY_v2_Count_QA = Form("v2_Corr_Count_Centrality_%d_%s",10,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_Count_QA] ,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.28,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.275,0.03,1,0.0,42,0);

  string KEY_v2_BW_QA = Form("v2_Corr_BW_Centrality_%d_%s",10,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_BW_QA] ,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.255,0.03,1,0.0,42,0);

  //string leg_besi = "BES-I";
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19_0,25,1,1.1);
  //Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,25,1,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.235,0.03,1,0.0,42,0);


  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.275,0.04,1,0.0,42,0);
  string leg_centrality = "0%-10%";
  plotTopLegend((char*)leg_centrality.c_str(),3.1,0.245,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/v2_pT_%s_0-10.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_v2_pT->SaveAs(figure_name.c_str());
  }
  {
  TCanvas *c_v2_pT = new TCanvas("c_v2_pT_2","c_v2_pT_2",10,10,800,800);
  c_v2_pT->cd();
  c_v2_pT->cd()->SetLeftMargin(0.15);
  c_v2_pT->cd()->SetBottomMargin(0.15);
  c_v2_pT->cd()->SetTicks(1,1);
  c_v2_pT->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(0.0,0.3);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");

  string KEY_v2_Count_QA = Form("v2_Corr_Count_Centrality_%d_%s",11,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_Count_QA] ,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.28,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.275,0.03,1,0.0,42,0);

  string KEY_v2_BW_QA = Form("v2_Corr_BW_Centrality_%d_%s",11,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_BW_QA] ,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.255,0.03,1,0.0,42,0);

  //string leg_besi = "BES-I";
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19_10,25,1,1.1);
  //Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,25,1,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.235,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.275,0.04,1,0.0,42,0);
  string leg_centrality = "10%-40%";
  plotTopLegend((char*)leg_centrality.c_str(),3.1,0.245,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/v2_pT_%s_10-40.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_v2_pT->SaveAs(figure_name.c_str());
  }
  {
  TCanvas *c_v2_pT = new TCanvas("c_v2_pT_3","c_v2_pT_3",10,10,800,800);
  c_v2_pT->cd();
  c_v2_pT->cd()->SetLeftMargin(0.15);
  c_v2_pT->cd()->SetBottomMargin(0.15);
  c_v2_pT->cd()->SetTicks(1,1);
  c_v2_pT->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(0.0,0.3);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");

  string KEY_v2_Count_QA = Form("v2_Corr_Count_Centrality_%d_%s",12,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_Count_QA] ,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.28,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.275,0.03,1,0.0,42,0);

  string KEY_v2_BW_QA = Form("v2_Corr_BW_Centrality_%d_%s",12,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2_corr[KEY_v2_BW_QA] ,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.255,0.03,1,0.0,42,0);

  //string leg_besi = "BES-I";
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19_40,25,1,1.1);
  //Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,25,1,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.235,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.275,0.04,1,0.0,42,0);
  string leg_centrality = "40%-80%";
  plotTopLegend((char*)leg_centrality.c_str(),3.1,0.245,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/v2_pT_%s_40-80.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_v2_pT->SaveAs(figure_name.c_str());
  }
  
  string OutPutFile = "../output/RawPhiv2/RawPhiV2PtSysKStar.root"; 
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
  {
    string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
    g_mv2[KEY_v2_Count]->Write();
    string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
    g_mv2[KEY_v2_BW]->Write();
    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      string KEY_Count = Form("v2_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
      h_mCounts[KEY_Count]->Write();
      string KEY_BW = Form("v2_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
      h_mCounts[KEY_BW]->Write();
    }
  }
  File_OutPut->Close();
}
