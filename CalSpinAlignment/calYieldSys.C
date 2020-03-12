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
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calYieldSys(int energy = 3, int pid = 0)
{
  string InPutFile_SE = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/Yields/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  // string InPutFile_SE = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  string InPutFile_ME = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/Yields/Yields_ME_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  // string InPutFile_ME = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/Yields_ME_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in pT spectra histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMassSpec_SE, h_mMassSpec_ME, h_mMassSpec_SM;
  string pT[2] = {"low","high"};
  for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(int i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(int i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
	  for(int i_range = 0; i_range < 2; i_range++)
	  {
	    string KEY_InPutSE = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_SE",i_pt,pT[i_range].c_str(),i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str());
	    h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 

	    string KEY_InPutME = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_ME",i_pt,pT[i_range].c_str(),i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str());
	    h_mInPut_ME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 

	    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	    {
	      string KEY_SE = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,pT[i_range].c_str(),i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMassSpec_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());

	      string KEY_ME = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_pt,pT[i_range].c_str(),i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMassSpec_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_ME.c_str());

	      string KEY_SM = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_pt,pT[i_range].c_str(),i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMassSpec_SM[KEY_SM] = (TH1F*)h_mMassSpec_SE[KEY_SE]->Clone();

	      if(i_norm < 2)
	      {
		int Norm_bin_start = h_mMassSpec_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_norm]);
		int Norm_bin_stop  = h_mMassSpec_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_norm]);

		float Inte_SE = h_mMassSpec_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
		float Inte_ME = h_mMassSpec_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);

		if(Inte_ME > 0.0)
		{
		  h_mMassSpec_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
		  h_mMassSpec_SM[KEY_SM]->Add(h_mMassSpec_ME[KEY_ME],-1.0);
		}
	      }
	      else
	      {
		float Inte_SE = 0.0;
		float Inte_ME = 0.0;

		for(int i_inte = 0; i_inte < 2; ++i_inte)
		{
		  int Norm_bin_start = h_mMassSpec_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_inte]);
		  int Norm_bin_stop  = h_mMassSpec_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_inte]);
		  Inte_SE += h_mMassSpec_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
		  Inte_ME += h_mMassSpec_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
		}

		if(Inte_ME > 0.0)
		{
		  h_mMassSpec_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
		  h_mMassSpec_SM[KEY_SM]->Add(h_mMassSpec_ME[KEY_ME],-1.0);
		}
	      }
	    }
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

  string KEY_SE_QA = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SE",vmsa::pt_RawQA[energy],pT[0].c_str(),vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMassSpec_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMassSpec_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMassSpec_SE[KEY_SE_QA]->GetMaximum());
  h_mMassSpec_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_ME",vmsa::pt_RawQA[energy],pT[0].c_str(),vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMassSpec_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMassSpec_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMassSpec_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMassSpec_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SM",vmsa::pt_RawQA[energy],pT[0].c_str(),vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMassSpec_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMassSpec_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMassSpec_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMassSpec_SM[KEY_SM_QA]->DrawCopy("h same");

  if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
  {
    PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMassSpec_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMassSpec_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
  }
  if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
  {
    PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMassSpec_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMassSpec_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
  }
#endif

#if _PlotQA_
  // QA Plots for pT bins
  TCanvas *c_pT_low = new TCanvas("c_pT_low","c_pT_low",10,10,1400,1400);
  c_pT_low->Divide(5,5);
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    c_pT_low->cd(i_pt+1);
    c_pT_low->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT_low->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT_low->cd(i_pt+1)->SetTicks(1,1);
    c_pT_low->cd(i_pt+1)->SetGrid(0,0);

    string KEY_SE_QA = Form("Spec_pt_%d_low_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec_SE[KEY_SE_QA]->DrawCopy("PE");

    string KEY_ME_QA = Form("Spec_pt_%d_low_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMassSpec_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMassSpec_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMassSpec_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("Spec_pt_%d_low_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMassSpec_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMassSpec_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMassSpec_SM[KEY_SM_QA]->DrawCopy("h same");

    TString pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],0.5*(vmsa::ptRawStart[i_pt]+vmsa::ptRawStop[i_pt]));
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
  }

  TCanvas *c_pT_high = new TCanvas("c_pT_high","c_pT_high",10,10,1400,1400);
  c_pT_high->Divide(5,5);
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    c_pT_high->cd(i_pt+1);
    c_pT_high->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT_high->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT_high->cd(i_pt+1)->SetTicks(1,1);
    c_pT_high->cd(i_pt+1)->SetGrid(0,0);

    string KEY_SE_QA = Form("Spec_pt_%d_high_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec_SE[KEY_SE_QA]->DrawCopy("PE");

    string KEY_ME_QA = Form("Spec_pt_%d_high_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMassSpec_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMassSpec_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMassSpec_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("Spec_pt_%d_high_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMassSpec_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMassSpec_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMassSpec_SM[KEY_SM_QA]->DrawCopy("h same");

    TString pT_range = Form("[%.2f,%.2f]",0.5*(vmsa::ptRawStart[i_pt]+vmsa::ptRawStop[i_pt]),vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  // pT rebin
  TH1FMap h_mMassSpec; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	  {
	    string KEY = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d",pt_bin,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	    // cout << endl;
	    // cout << "KEY = " << KEY.c_str() << " => " << endl;
	    for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	    {
	      string KEY_SM_low = Form("Spec_pt_%d_low_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      string KEY_SM_high = Form("Spec_pt_%d_high_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      // cout << KEY_SM_low.c_str() << " + " << KEY_SM_high.c_str() << endl;
	      if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
	      {
		h_mMassSpec[KEY] = (TH1F*)h_mMassSpec_SM[KEY_SM_low]->Clone(KEY.c_str());
		h_mMassSpec[KEY]->Add(h_mMassSpec_SM[KEY_SM_high],1.0);
	      }
	      else
	      {
		h_mMassSpec[KEY]->Add(h_mMassSpec_SM[KEY_SM_low],1.0);
		h_mMassSpec[KEY]->Add(h_mMassSpec_SM[KEY_SM_high],1.0);
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
    string KEY_QA = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec[KEY_QA]->SetMarkerStyle(24);
    h_mMassSpec[KEY_QA]->SetMarkerSize(1.2);
    h_mMassSpec[KEY_QA]->SetLineColor(kGray+2);
    h_mMassSpec[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }
#endif

  // subtract linear background
  // Poly + Breit Wignar fit to InvMass
  vecFMap ParFit;
  TH1FMap h_mMassSpec_QA;

  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(int i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(int i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    string KEY = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMassSpec_QA[KEY] = (TH1F*)h_mMassSpec[KEY]->Clone();
	    TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
	    for(int i_par = 0; i_par < 5; i_par++)
	    {
	      f_bw->ReleaseParameter(i_par);
	    }
	    f_bw->SetParameter(0,vmsa::InvMass[pid]);
	    f_bw->SetParLimits(0,vmsa::InvMass[pid]-1.5*vmsa::Width[pid],vmsa::InvMass[pid]+1.5*vmsa::Width[pid]);
	    f_bw->SetParameter(1,vmsa::Width[pid]);
	    f_bw->SetParLimits(1,0.004,0.070);
	    f_bw->SetParameter(2,1.0);
	    f_bw->SetParameter(3,-1.0);
	    f_bw->SetParameter(4,1.0);
	    f_bw->SetParameter(2,h_mMassSpec[KEY]->GetMaximum()/f_bw->GetMaximum());
	    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	    ParFit[KEY].clear();
	    // cout << "i_pt = " << i_pt << ", i_cent = " << i_cent << ", i_dca = " << i_dca << ", i_sig = " << i_sig << ", i_norm = " << i_norm << endl;
	    h_mMassSpec[KEY]->Fit(f_bw,"MQNR");
	    for(int n_par = 0; n_par < 5; n_par++)
	    {
	      ParFit[KEY].push_back(static_cast<float>(f_bw->GetParameter(n_par)));
	    }

	    TF1 *f_poly = new TF1("f_poly",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2); // poly fit for linear background
	    f_poly->SetParameter(0,ParFit[KEY][3]);
	    f_poly->SetParameter(1,ParFit[KEY][4]);

	    h_mMassSpec[KEY]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
  TCanvas *c_mMassSpec = new TCanvas("c_mMassSpec","c_mMassSpec",10,10,1400,1400);
  c_mMassSpec->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_mMassSpec->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_mMassSpec->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mMassSpec->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mMassSpec->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_mMassSpec->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);

    string KEY_QA = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec[KEY_QA]->SetMarkerColor(1);
    h_mMassSpec[KEY_QA]->SetMarkerStyle(24);
    h_mMassSpec[KEY_QA]->SetMarkerSize(0.8);
    h_mMassSpec[KEY_QA]->DrawCopy("PE");

    h_mMassSpec_QA[KEY_QA]->SetMarkerColor(4);
    h_mMassSpec_QA[KEY_QA]->SetMarkerStyle(24);
    h_mMassSpec_QA[KEY_QA]->SetMarkerSize(0.8);
    h_mMassSpec_QA[KEY_QA]->DrawCopy("PE same");

    TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
    for(int i_par = 0; i_par < 5; i_par++)
    {
      f_bw->SetParameter(i_par,ParFit[KEY_QA][i_par]);
    }
    f_bw->SetLineColor(2);
    f_bw->SetLineStyle(1);
    f_bw->SetLineWidth(2);
    f_bw->DrawCopy("l same");

    TF1 *f_poly = new TF1("f_poly",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
    f_poly->SetParameter(0,ParFit[KEY_QA][3]);
    f_poly->SetParameter(1,ParFit[KEY_QA][4]);
    f_poly->SetLineColor(2);
    f_poly->SetLineStyle(2);
    f_poly->SetLineWidth(4);
    f_poly->DrawCopy("l same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }
#endif

  // calculate total yields for each pT bin via gaussian and breit wigner fits
  // if statistics are not enough => set yields to -1.0+/-0.5
  vecFMap ParYield;
  vecFMap yields_count, yields_inte;
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(int i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(int i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    string KEY = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
	    f_yields_bw->SetParameter(0,vmsa::InvMass[pid]);
	    f_yields_bw->SetParLimits(0,vmsa::InvMass[pid]-0.005,vmsa::InvMass[pid]+0.005);
	    f_yields_bw->SetParameter(1,vmsa::Width[pid]);
	    f_yields_bw->SetParameter(2,h_mMassSpec[KEY]->GetMaximum());
	    f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);

	    int NumofBin = h_mMassSpec[KEY]->FindBin(vmsa::InvMass[pid]);
	    if(h_mMassSpec[KEY]->GetBinContent(NumofBin) > 0.0)
	    {
	      h_mMassSpec[KEY]->Fit(f_yields_bw,"MQNR");
	      ParYield[KEY].clear();
	      ParYield[KEY].push_back(static_cast<float>(f_yields_bw->GetParameter(0)));
	      ParYield[KEY].push_back(static_cast<float>(f_yields_bw->GetParameter(1)));
	      ParYield[KEY].push_back(static_cast<float>(f_yields_bw->GetParameter(2)));

	      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	      {
		// bin counting
		float counts_gaus = 0.0;
		float errors_gaus = 0.0;
		int bin_start = h_mMassSpec[KEY]->FindBin(ParYield[KEY][0]-vmsa::nSigVecSys[i_sigma]*ParYield[KEY][1]);
		int bin_stop  = h_mMassSpec[KEY]->FindBin(ParYield[KEY][0]+vmsa::nSigVecSys[i_sigma]*ParYield[KEY][1]);
		for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
		{
		  counts_gaus += h_mMassSpec[KEY]->GetBinContent(i_bin);
		  errors_gaus += h_mMassSpec[KEY]->GetBinError(i_bin)*h_mMassSpec[KEY]->GetBinError(i_bin);
		}
		string KEY_yields = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
		yields_count[KEY_yields].clear();
		yields_count[KEY_yields].push_back(static_cast<float>(counts_gaus));
		yields_count[KEY_yields].push_back(static_cast<float>(TMath::Sqrt(errors_gaus)));

		// integrating for breit wigner
		float bin_width = h_mMassSpec[KEY]->GetBinWidth(1);
		float Inte_start = ParYield[KEY][0]-vmsa::nSigVecSys[i_sigma]*ParYield[KEY][1]-0.5*bin_width;
		float Inte_stop  = ParYield[KEY][0]+vmsa::nSigVecSys[i_sigma]*ParYield[KEY][1]+0.5*bin_width;
		float counts_bw = f_yields_bw->Integral(Inte_start,Inte_stop)/bin_width;
		float errors_bw = f_yields_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
		yields_inte[KEY_yields].clear();
		yields_inte[KEY_yields].push_back(static_cast<float>(counts_bw));
		yields_inte[KEY_yields].push_back(static_cast<float>(errors_bw));
	      }
	    }
	    else
	    {
	      ParYield[KEY].clear();
	      ParYield[KEY].push_back(static_cast<float>(f_yields_bw->GetParameter(0)));
	      ParYield[KEY].push_back(static_cast<float>(f_yields_bw->GetParameter(1)));
	      ParYield[KEY].push_back(static_cast<float>(f_yields_bw->GetParameter(2)));

	      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	      {
		string KEY_yields = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
		yields_count[KEY_yields].clear();
		yields_count[KEY_yields].push_back(static_cast<float>(-1.0));
		yields_count[KEY_yields].push_back(static_cast<float>(0.5));

		yields_inte[KEY_yields].clear();
		yields_inte[KEY_yields].push_back(static_cast<float>(-1.0));
		yields_inte[KEY_yields].push_back(static_cast<float>(0.5));
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA: different counting method: bin counting vs breit wigner integrating
  TCanvas *c_Yields = new TCanvas("c_Yields","c_Yields",10,10,1400,1400);
  c_Yields->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_Yields->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_Yields->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_Yields->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_Yields->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_Yields->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);

    string KEY_QA = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMassSpec[KEY_QA]->SetTitle("");
    h_mMassSpec[KEY_QA]->SetStats(0);
    h_mMassSpec[KEY_QA]->SetMarkerStyle(24);
    h_mMassSpec[KEY_QA]->SetMarkerColor(kGray+3);
    h_mMassSpec[KEY_QA]->SetMarkerSize(0.8);
    h_mMassSpec[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mMassSpec[KEY_QA]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mMassSpec[KEY_QA]->GetXaxis()->CenterTitle();
    h_mMassSpec[KEY_QA]->GetYaxis()->SetTitle("Counts");
    h_mMassSpec[KEY_QA]->GetYaxis()->CenterTitle();
    h_mMassSpec[KEY_QA]->GetYaxis()->SetTitleOffset(1.2);
    h_mMassSpec[KEY_QA]->DrawCopy("pE");

    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
    f_yields_bw->SetParameter(0,ParYield[KEY_QA][0]);
    f_yields_bw->SetParameter(1,ParYield[KEY_QA][1]);
    f_yields_bw->SetParameter(2,ParYield[KEY_QA][2]);
    f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    f_yields_bw->SetLineColor(2);
    f_yields_bw->SetLineStyle(1);
    f_yields_bw->SetLineWidth(2);
    f_yields_bw->DrawCopy("l same");

    float x1 = ParYield[KEY_QA][0] - vmsa::nSigVecSys[vmsa::Sig_QA]*ParYield[KEY_QA][1];
    float x2 = ParYield[KEY_QA][0] + vmsa::nSigVecSys[vmsa::Sig_QA]*ParYield[KEY_QA][1];
    float y = h_mMassSpec[KEY_QA]->GetBinContent(h_mMassSpec[KEY_QA]->FindBin(ParYield[KEY_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }
#endif

  // declare histogram with different pT width
  int const NumOfBin = vmsa::pt_rebin_last[energy]+1;
  float ptBin_low[NumOfBin], ptBin_up[NumOfBin], ptBin_width[NumOfBin];
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    ptBin_low[i_pt] = vmsa::pt_low[energy][i_pt];
    ptBin_up[i_pt]  = vmsa::pt_up[energy][i_pt];
    ptBin_width[i_pt] = vmsa::pt_up[energy][i_pt]-vmsa::pt_low[energy][i_pt];
    // cout << "pt_start = " << ptBin_low[i_pt] << ", pt_stop = " << ptBin_up[i_pt] << ", width = " << ptBin_width[i_pt] << endl;
  }
  ptBin_low[NumOfBin-1] = vmsa::pt_low[energy][NumOfBin-1];
  ptBin_up[NumOfBin-1] = vmsa::pt_up[energy][NumOfBin-1];
  ptBin_width[NumOfBin-1] = vmsa::pt_up[energy][NumOfBin-1]-vmsa::pt_low[energy][NumOfBin-1];
  // cout << "pt_start = " << ptBin_low[NumOfBin-1] << ", pt_stop = " << ptBin_up[NumOfBin-1] << ", width = " << ptBin_width[NumOfBin-1] << endl;

  TH1FMap h_mPt;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  {
    for(int i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(int i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    string KEY_pT_counts = Form("Yield_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_Count",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
	    h_mPt[KEY_pT_counts] = new TH1F(KEY_pT_counts.c_str(),KEY_pT_counts.c_str(),vmsa::pt_rebin_last[energy],ptBin_low);
	    string KEY_pT_inte = Form("Yield_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_Inte",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
	    h_mPt[KEY_pT_inte] = new TH1F(KEY_pT_inte.c_str(),KEY_pT_inte.c_str(),vmsa::pt_rebin_last[energy],ptBin_low);
	    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
	    {
	      string KEY_yields = Form("Yield_pt_%d_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
	      h_mPt[KEY_pT_counts]->SetBinContent(i_pt+1,yields_count[KEY_yields][0]/ptBin_width[i_pt]);
	      h_mPt[KEY_pT_counts]->SetBinError(i_pt+1,yields_count[KEY_yields][1]/TMath::Sqrt(ptBin_width[i_pt]));
	      h_mPt[KEY_pT_inte]->SetBinContent(i_pt+1,yields_inte[KEY_yields][0]/ptBin_width[i_pt]);
	      h_mPt[KEY_pT_inte]->SetBinError(i_pt+1,yields_inte[KEY_yields][1]/TMath::Sqrt(ptBin_width[i_pt]));
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA: pt spectra
  TCanvas *c_Pt = new TCanvas("c_Pt","c_Pt",10,10,800,800);
  c_Pt->cd();
  c_Pt->cd()->SetLeftMargin(0.20);
  c_Pt->cd()->SetBottomMargin(0.20);
  c_Pt->cd()->SetTicks(1,1);
  c_Pt->cd()->SetGrid(0,0);
  // c_Pt->cd()->SetLogy();
  string KEY_pT_counts_QA = Form("Yield_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_Count",vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Sig_QA);
  TH1F *h_play = new TH1F("h_play","h_play",110,-1.0,10.0);
  for(Int_t i_bin = 0; i_bin < 110; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-1000.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitle("dN/dp_{T}");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitleOffset(1.2);
  h_play->GetYaxis()->SetRangeUser(-5.0,1.2*h_mPt[KEY_pT_counts_QA]->GetMaximum());
  h_play->DrawCopy("pE");
  h_mPt[KEY_pT_counts_QA]->SetMarkerStyle(24);
  h_mPt[KEY_pT_counts_QA]->SetMarkerColor(4);
  h_mPt[KEY_pT_counts_QA]->SetMarkerSize(1.0);
  h_mPt[KEY_pT_counts_QA]->DrawCopy("pE same");

  string KEY_pT_inte_QA = Form("Yield_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_Inte",vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Sig_QA);
  h_mPt[KEY_pT_inte_QA]->SetMarkerColor(2);
  h_mPt[KEY_pT_inte_QA]->SetMarkerStyle(24);
  h_mPt[KEY_pT_inte_QA]->SetMarkerSize(1.0);
  h_mPt[KEY_pT_inte_QA]->DrawCopy("pE same");
#endif

  string outputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/RawYieldSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawYieldSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  {
    for(int i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(int i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    string KEY_pT_counts = Form("Yield_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_Count",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
	    h_mPt[KEY_pT_counts]->Write();
	    string KEY_pT_inte = Form("Yield_Centrality_%d_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_Inte",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma);
	    h_mPt[KEY_pT_inte]->Write();
	  }
	}
      }
    }
  }
  File_OutPut->Close();

}
