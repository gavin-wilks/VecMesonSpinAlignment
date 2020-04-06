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
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

#ifndef _SaveQA_
#define _SaveQA_  0
#endif

void calSpinAlignmentSys_SideBand(int energy = 6, int pid = 0)
{
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  string InPutFile_SE = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/Yields/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  // string InPutFile_SE = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  int norm = vmsa::Norm_QA;
  // read in histogram for same event
  TH1FMap h_mInPut_SE;
  TH1FMap h_mMass_SE;
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
	    string KEY_InPutSE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_SE",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str());
	    h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 

	    string KEY_SE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	    h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
  c_peak->cd();
  c_peak->cd()->SetLeftMargin(0.15);
  c_peak->cd()->SetBottomMargin(0.15);
  c_peak->cd()->SetTicks(1,1);
  c_peak->cd()->SetGrid(0,0);
  string KEY_SE_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_SE",vmsa::pt_RawQA[energy],vmsa::Cent_start,vmsa::CTS_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),norm);
  h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  PlotLine(vmsa::Norm_Start[pid][norm],vmsa::Norm_Start[pid][norm],0,0.8*h_mMass_SE[KEY_SE_QA]->GetMaximum(),4,2,2);
  PlotLine(vmsa::Norm_Stop[pid][norm],vmsa::Norm_Stop[pid][norm],0,0.8*h_mMass_SE[KEY_SE_QA]->GetMaximum(),4,2,2);


  // QA Plots for pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
  c_pT->Divide(5,5);
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    c_pT->cd(i_pt+1);
    c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT->cd(i_pt+1)->SetTicks(1,1);
    c_pT->cd(i_pt+1)->SetGrid(0,0);
    string KEY_SE_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,vmsa::Cent_start,vmsa::CTS_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),norm);
    h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
    h_mMass_SE[KEY_SE_QA]->DrawCopy();

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);

    PlotLine(vmsa::Norm_Start[pid][norm],vmsa::Norm_Start[pid][norm],0,0.8*h_mMass_SE[KEY_SE_QA]->GetMaximum(),4,2,2);
    PlotLine(vmsa::Norm_Stop[pid][norm],vmsa::Norm_Stop[pid][norm],0,0.8*h_mMass_SE[KEY_SE_QA]->GetMaximum(),4,2,2);
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
	  for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",pt_bin,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	    for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	    {
	      string KEY_SE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	      // cout << "KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	      if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
	      {
		h_mMass[KEY] = (TH1F*)h_mMass_SE[KEY_SE]->Clone(KEY.c_str());
	      }
	      else
	      {
		h_mMass[KEY]->Add(h_mMass_SE[KEY_SE],1.0);
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for pT rebins
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::CTS_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),norm);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(1.2);
    h_mMass[KEY_QA]->SetLineColor(kGray+2);
    h_mMass[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }
#endif

  // extract counts vs. pT with diffenretial integration ranges and methods
  TH1FMap h_mCounts;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
	string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	g_mRho[KEY_rho] = new TGraphAsymmErrors();
	for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	{
	  string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	  h_mCounts[KEY_counts] = new TH1F(KEY_counts.c_str(),KEY_counts.c_str(),7,0.0,1.0);

	  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; ++i_theta) // cos(theta*) loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	    float bin_center = 1/14.0+i_theta/7.0;
	    int bin_start = h_mMass[KEY]->FindBin(vmsa::Norm_Start[pid][norm]);
	    int bin_stop  = h_mMass[KEY]->FindBin(vmsa::Norm_Stop[pid][norm]);
	    float counts = 0.0;
	    float errors = 0.0;
	    for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	    {
	      counts += h_mMass[KEY]->GetBinContent(i_bin);
	      errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	    }
	    h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts);
	    h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors));
	  }
	  float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

	  TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
	  f_rho->SetParameter(0,0.33);
	  f_rho->SetParameter(1,h_mCounts[KEY_counts]->GetMaximum());
	  h_mCounts[KEY_counts]->Fit(f_rho,"NMRI");
	  Par_rhoFit[KEY_counts].clear();
	  Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
	  Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
	  g_mRho[KEY_rho]->SetPoint(i_pt,pt_mean,f_rho->GetParameter(0));
	  g_mRho[KEY_rho]->SetPointError(i_pt,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
	}
      }
    }
  }
  
#if _PlotQA_
  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,900);
  c_diff->Divide(3,3);
  for(int i_theta = 0; i_theta < 9; ++i_theta)
  {
    c_diff->cd(i_theta+1);
    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
    c_diff->cd(i_theta+1)->SetTicks(1,1);
    c_diff->cd(i_theta+1)->SetGrid(0,0);
    if(i_theta < vmsa::CTS_stop)
    {
      string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",vmsa::pt_QA[energy],vmsa::Cent_start,i_theta,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),norm);
      h_mMass[KEY_QA]->SetTitle("");
      h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass[KEY_QA]->GetXaxis()->CenterTitle();

      h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
      h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
      h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
      h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
      h_mMass[KEY_QA]->GetYaxis()->CenterTitle();

      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
      h_mMass[KEY_QA]->SetMarkerSize(1.2);
      h_mMass[KEY_QA]->Draw("pE");
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
      PlotLine(vmsa::Norm_Start[pid][norm],vmsa::Norm_Start[pid][norm],0,0.95*h_mMass[KEY_QA]->GetMaximum(),4,2,2);
      PlotLine(vmsa::Norm_Stop[pid][norm],vmsa::Norm_Stop[pid][norm],0,0.95*h_mMass[KEY_QA]->GetMaximum(),4,2,2);
    }
  }

  c_diff->cd(9);
  string KEY_counts_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",vmsa::pt_QA[energy],vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),norm);
  h_mCounts[KEY_counts_QA]->SetTitle("");
  h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
  h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
  h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta*)");
  h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
  h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();

  h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
  h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
  h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("Counts");
  h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
  h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
  h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();

  h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
  h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
  h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
  h_mCounts[KEY_counts_QA]->Draw("pE");

  TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
  f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
  f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
  f_rho->SetLineWidth(2);
  f_rho->SetLineStyle(2);
  f_rho->Draw("l same");
  // c_diff->SaveAs("../figures/c_diff_2.eps");
#endif

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.2,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
	string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_dca+10*i_sig+1,1.1);
      }
    }
  }
  // c_rho->SaveAs("../figures/c_rho.eps");

  string outputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/RawRhoPtSys_SideBand.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
	string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
	g_mRho[KEY_rho]->Write();
	for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	{
	  string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),norm);
	  h_mCounts[KEY_counts]->Write();
	}
      }
    }
  }
  File_OutPut->Close();
}
