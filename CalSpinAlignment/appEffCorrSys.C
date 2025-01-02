#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

using namespace std;

int const Color_Counts = kGray+2;
int const Color_BW     = kAzure-4;

TH1F* CalEffCorr(TH1F *h_counts, TH1D *h_eff, std::string HistName)
{
  gStyle->SetOptDate(0);
  TH1F* h_effCorr = (TH1F*)h_counts->Clone();
  h_effCorr->Divide(h_eff);
  for(int i_bin = 1; i_bin < h_effCorr->GetNbinsX()+1; ++i_bin)
  {
    float counts = h_counts->GetBinContent(i_bin);
    float err_counts = h_counts->GetBinError(i_bin);
    float eff = h_eff->GetBinContent(i_bin);
    float err_eff = h_eff->GetBinError(i_bin);
    float err = ErrDiv(counts,eff,err_counts,err_eff);
    h_effCorr->SetBinError(i_bin,err);
  }
  h_effCorr->SetName(HistName.c_str());

  return h_effCorr;
}

TH1F* CalEffCorr(TH1F *h_counts, TF1 *f_eff, std::string HistName)
{
  TH1F* h_effCorr = (TH1F*)h_counts->Clone();
  h_effCorr->Divide(f_eff);
  float eff = f_eff->GetParameter(0);
  float err_eff = f_eff->GetParError(0);
  for(int i_bin = 1; i_bin < h_effCorr->GetNbinsX()+1; ++i_bin)
  {
    float counts = h_counts->GetBinContent(i_bin);
    float err_counts = h_counts->GetBinError(i_bin);
    float err = ErrDiv(counts,eff,err_counts,err_eff);
    h_effCorr->SetBinError(i_bin,err);
  }
  h_effCorr->SetName(HistName.c_str());
  h_effCorr->SetTitle(HistName.c_str());

  return h_effCorr;
}

// used for 20%-60% and eta_gap = 0.05 only
void appEffCorrSys(int energy = 4, int pid = 0, int ptQA = 2, bool isBesI = false, int order = 2, bool random3D = false, int etamode = 5)
{
  std::string etastring;
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 1) etastring = "eta1_eta1p5";
  if(etamode == 2) etastring = "eta1p5_eta1p5";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";

  string InPutRho = Form("../output/AuAu%s/%s/RawPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str());
  if(order == 1) InPutRho = Form("../output/AuAu%s/%s/RawPhiPtSys_Psi1.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(random3D) InPutRho = Form("../output/AuAu%s/%s/3DRandom/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(isBesI) InPutRho = Form("../data/BESI/RawRhoPtSys.root");
  TFile *File_Rho = TFile::Open(InPutRho.c_str());
  cout << InPutRho << endl;

  string InPutEff = Form("./efficiencies/%s/%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode);
  if(isBesI) InPutEff = Form("../../../FileTransfers/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str());
  cout << InPutEff << endl;

  TFile *File_Eff = TFile::Open(InPutEff.c_str());

  TH1FMap h_mCounts, h_mCountsEff;
  TH1DMap h_mEff;
  TGraMap g_mRho;
  vecFMap Par_rhoFit;
  for(int i_cent = 9;/*vmsa::Cent_start;*/ i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              cout << KEY_rho << endl;
	      g_mRho[KEY_rho] = new TGraphAsymmErrors();
	      for(int i_pt = /*vmsa::pt_rebin_first[energy]*/2; i_pt < 6/*vmsa::pt_rebin_last[energy]*/; ++i_pt) // pt loop
	      {
		string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
                cout << KEY_counts << endl;
		h_mCounts[KEY_counts] = (TH1F*)File_Rho->Get(KEY_counts.c_str());
		h_mCounts[KEY_counts]->Print();
		string KEY_Eff = Form("h_mEffCos_Cent_9_Pt_%d",i_pt); // 20%-60%
                cout << KEY_Eff << endl;
		h_mEff[KEY_Eff] = (TH1D*)File_Eff->Get(KEY_Eff.c_str());
		TF1 *f_poly = new TF1("f_poly","pol0",0.0,1.0);
		h_mEff[KEY_Eff]->Fit(f_poly,"NQ");
		h_mEff[KEY_Eff]->Print();
                
                //cout << "Did I fit?" << endl;
	        h_mCountsEff[KEY_counts] = CalEffCorr(h_mCounts[KEY_counts],h_mEff[KEY_Eff],KEY_counts);
                //cout << "Calculated the efficiency corrected value" << endl;
  		//h_mCountsEff[KEY_counts] = CalEffCorr(h_mCounts[KEY_counts],f_poly,KEY_counts);
    
		float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

		TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
		f_rho->SetParameter(0,0.33);
		f_rho->SetParameter(1,h_mCountsEff[KEY_counts]->GetMaximum());
		h_mCountsEff[KEY_counts]->Fit(f_rho,"NMRI");
		Par_rhoFit[KEY_counts].clear();
		Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
		Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
		g_mRho[KEY_rho]->SetPoint(i_pt,pt_mean,f_rho->GetParameter(0));
		g_mRho[KEY_rho]->SetPointError(i_pt,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
	      }
	    }
	  }
	}
      }
    }
  }

/*#if _PlotQA_

    double binpt[7] = {0.4,0.8,1.2,1.8,2.4,3.0,4.2};
    TCanvas *c_Counts = new TCanvas("c_Counts","c_Counts",10,10,1000,800);
    c_Counts->Divide(2,2);
                
    for(int ipt = 2; ipt < 6; ipt++)
    { 
      c_Counts->cd(ipt-1);
      c_Counts->cd(ipt-1)->SetLeftMargin(0.15);
      c_Counts->cd(ipt-1)->SetBottomMargin(0.15);
      c_Counts->cd(ipt-1)->SetTicks(1,1);
      c_Counts->cd(ipt-1)->SetGrid(0,0);

      string KEY_counts_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ipt,9,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      h_mCountsEff[KEY_counts_QA];
      h_mCountsEff[KEY_counts_QA]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c",binpt[ipt],binpt[ipt+1]));
      h_mCountsEff[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mCountsEff[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mCountsEff[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta*)");
      h_mCountsEff[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mCountsEff[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mCountsEff[KEY_counts_QA]->GetXaxis()->CenterTitle();
      h_mCountsEff[KEY_counts_QA]->SetStats(0);

      h_mCountsEff[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCountsEff[KEY_counts_QA]->GetMaximum());
      h_mCountsEff[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
      h_mCountsEff[KEY_counts_QA]->GetYaxis()->SetTitle("Counts");
      h_mCountsEff[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
      h_mCountsEff[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
      h_mCountsEff[KEY_counts_QA]->GetYaxis()->CenterTitle();

      h_mCountsEff[KEY_counts_QA]->SetMarkerStyle(20);
      h_mCountsEff[KEY_counts_QA]->SetMarkerColor(kRed);
      h_mCountsEff[KEY_counts_QA]->SetLineColor(kRed);
      h_mCountsEff[KEY_counts_QA]->SetLineWidth(1);
      h_mCountsEff[KEY_counts_QA]->SetMarkerSize(1.5);
      h_mCountsEff[KEY_counts_QA]->Draw("pE");

      TF1 *f_rhoE = new TF1("f_rhoE",SpinDensity,0.0,1.0,2);
      f_rhoE->SetParameter(0,0.33);
      f_rhoE->SetParameter(1,h_mCountsEff[KEY_counts_QA]->GetMaximum());
      h_mCountsEff[KEY_counts_QA]->Fit(f_rhoE,"NMRI");
      f_rhoE->SetLineColor(kBlue);
      f_rhoE->Draw("same");

      h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
      h_mCounts[KEY_counts_QA]->SetTitle(Form("%.1f < p_{T} < %.1f GeV/c",binpt[ipt],binpt[ipt+1]));
      h_mCounts[KEY_counts_QA]->SetMarkerColor(kGray+1);
      h_mCounts[KEY_counts_QA]->SetLineColor(kGray+1);
      h_mCounts[KEY_counts_QA]->SetLineWidth(1);
      h_mCounts[KEY_counts_QA]->SetMarkerSize(1.5);
      h_mCounts[KEY_counts_QA]->SetStats(0);
      h_mCounts[KEY_counts_QA]->Draw("pE same");
     
      TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
      f_rho->SetParameter(0,0.33);
      f_rho->SetParameter(1,h_mCounts[KEY_counts_QA]->GetMaximum());
      h_mCounts[KEY_counts_QA]->Fit(f_rho,"NMRI");
      f_rho->SetLineColor(kBlack);
      f_rho->SetLineStyle(2);
      f_rho->Draw("same");

      TLegend *raw = new TLegend(0.2,0.25,0.8,0.4);
      string leg_raw = Form("#rho_{00}^{raw} = %.4f #pm %.4f",f_rho->GetParameter(0),f_rho->GetParError(0));
      raw->AddEntry(f_rho,leg_raw.c_str());
      raw->SetTextColor(kGray+1);
      raw->SetBorderSize(0);
      raw->SetFillColor(0);
      raw->Draw("same");      
     
      TLegend *eff = new TLegend(0.2,0.5,0.8,0.65);
      string leg_eff = Form("#rho_{00}^{obs} = %.4f #pm %.4f",f_rhoE->GetParameter(0),f_rhoE->GetParError(0));
      eff->AddEntry(f_rhoE,leg_eff.c_str());
      eff->SetTextColor(kRed);
      eff->SetBorderSize(0);
      eff->SetFillColor(0);     
      eff->Draw("same");

      c_Counts->SaveAs("figures/Preliminary/yieldsvscos.pdf");
    }
#endif
*/
  /*TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
*/  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
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
 /* h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i_cent = 9; i_cent < 10; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
	    }
	  }
	}
      }
    }
  }*/

  string OutPutFile = Form("../output/AuAu%s/%s/EffPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str());
  if(random3D) OutPutFile = Form("../output/AuAu%s/%s/3DRandom/EffPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(isBesI) OutPutFile = Form("../output/AuAu%s/%s/BESI/EffPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();

  for(int i_cent = 9; /*vmsa::Cent_start;*/ i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
	      g_mRho[KEY_rho]->Write();
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}
