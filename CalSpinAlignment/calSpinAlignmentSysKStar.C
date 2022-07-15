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

int Sig_start = 0;
int Sig_stop = 3;
float nSigVecSys[3] = {2.0,1.5,2.5};

void calSpinAlignmentSysKStar(int energy = 4, int pid = 2, int year = 0, bool random3D = false)
{
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  for(int i_pt = vmsa::pt_rebin_firstKS[energy]; i_pt < vmsa::pt_rebin_lastKS[energy]; i_pt++) // pt loop
  {
    for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < 10/*vmsa::Cent_stop*/; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
          for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++) // systematic loop for nSigma
          {
            if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	    for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	    {
	      string KEY_InteTheta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	      {
	        string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str());

	        if(i_theta == vmsa::CTS_start) h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)h_mMass[KEY]->Clone(KEY_InteTheta.c_str());
	        else h_mMass_InteTheta[KEY_InteTheta]->Add(h_mMass[KEY],1.0);
	      }
              cout << KEY_InteTheta << endl;
	      TF1 *f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
	      f_bw->SetParameter(0,vmsa::InvMass[pid]);
	      f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.01,vmsa::InvMass[pid]+0.01);
	      f_bw->SetParameter(1,vmsa::Width[pid]);
	      //f_bw->SetParLimits(1,0.01,0.2);
	      f_bw->SetParameter(2,1000.0);
	      float norm = h_mMass_InteTheta[KEY_InteTheta]->GetMaximum()/f_bw->GetMaximum();
	      f_bw->SetParameter(2,norm);
	      f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	      h_mMass_InteTheta[KEY_InteTheta]->Fit(f_bw,"MQNR");
	      Par_InteTheta[KEY_InteTheta].clear();
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(0)));
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(1)));
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(2)));
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(3)));
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(4)));
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(5)));
	    }
          }
	}
      }
    }
  }

/*#if _PlotQA_
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
      string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_0_%s_Norm_%d",vmsa::pt_QAKS[energy],9,i_theta,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
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
    }
    if(i_theta == vmsa::CTS_stop)
    {
      string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_0_%s_Norm_%d",vmsa::pt_QAKS[energy],9,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();

      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Yields");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();

      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
      TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
      f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][0]);
      f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][1]);
      f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][2]);
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(2);
      f_bw->SetLineWidth(2);
      f_bw->Draw("l same");
    }
  }
#endif*/

  /*
  // fit theta-differential bin to extract fit parameters for integration
  for(int i_pt = vmsa::pt_rebin_firstKS[energy]; i_pt < vmsa::pt_rebin_lastKS[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	{
	  for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	  {
	    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	    {
	      string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	    }
	  }
	}
      }
    }
  }
  */

  // extract counts vs. pT with diffenretial integration ranges and methods
  vecFMap Par;
  TH1FMap h_mCounts;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;
  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < 10/*vmsa::Cent_stop*/; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++) // systematic loop for nSigma
        {
          if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	  {
	    for(int i_sigma = Sig_start; i_sigma < Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
	        string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        g_mRho[KEY_rho] = new TGraphAsymmErrors();
	        for(int i_pt = vmsa::pt_rebin_firstKS[energy]; i_pt < vmsa::pt_rebin_lastKS[energy]; ++i_pt) // pt loop
	        {
	  	string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	  	h_mCounts[KEY_counts] = new TH1F(KEY_counts.c_str(),KEY_counts.c_str(),7,0.0,1.0);

	  	string KEY_InteTheta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	  	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; ++i_theta) // cos(theta*) loop
	  	{
	  	  string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
                  cout << KEY << endl;
	  	  TF1 *f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
	  	  f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta][0]);
	  	  f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta][1]);
	  	  f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta][2]/7.0);
	  	  f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta][3]/7.0);
	  	  f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta][4]/7.0);
	  	  f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta][5]);
	  	  f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  	  TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"MNRS");
                  
                    TF1 *f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
                    f_bg->SetParameter(0,f_bw->GetParameter(3));
                    f_bg->SetParameter(1,f_bw->GetParameter(4));
                    f_bg->SetParameter(2,f_bw->GetParameter(5));
                    f_bg->SetParError(0,f_bw->GetParError(3));
                    f_bg->SetParError(1,f_bw->GetParError(4));
                    f_bg->SetParError(2,f_bw->GetParError(5));

                    double params[3] = {result->GetParams()[3],result->GetParams()[4],result->GetParams()[5]};
                    TMatrixDSym covArr(3);
                    covArr(0,0) = result->GetCovarianceMatrix()(3,3);
                    covArr(0,1) = result->GetCovarianceMatrix()(3,4);
                    covArr(0,2) = result->GetCovarianceMatrix()(3,5);
                    covArr(1,0) = result->GetCovarianceMatrix()(4,3);
                    covArr(1,1) = result->GetCovarianceMatrix()(4,4);
                    covArr(1,2) = result->GetCovarianceMatrix()(4,5);
                    covArr(2,0) = result->GetCovarianceMatrix()(5,3);
                    covArr(2,1) = result->GetCovarianceMatrix()(5,4);
                    covArr(2,2) = result->GetCovarianceMatrix()(5,5);

	  	  float bin_width = h_mMass[KEY]->GetBinWidth(1);
	  	  float Inte_start = Par_InteTheta[KEY_InteTheta][0]-nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
                    float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
                    float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
                    float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width;

	  	  float bin_center = 1/14.0+i_theta/7.0;
	  	  if(i_method == 0)
	  	  {
	  	    int bin_start = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]-nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
	  	    int bin_stop  = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]+nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
	  	    float counts = 0.0;
	  	    float errors = 0.0;
	  	    for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	  	    {
	  	      counts += h_mMass[KEY]->GetBinContent(i_bin);
	  	      errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	  	    }
	  	    h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts-counts_bg);
	  	    h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors+errors_bg*errors_bg));
	  	  }
	  	  if(i_method == 1)
	  	  {
	  	    //float bin_width = h_mMass[KEY]->GetBinWidth(1);
	  	    //float Inte_start = Par_InteTheta[KEY_InteTheta][0]-nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
	  	    //float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
	  	    float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
	  	    float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	  	    h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts_bw-counts_bg);
	  	    h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));
	  	  }
	  	  Par[KEY].clear();
	  	  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(0)));
	  	  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(1)));
	  	  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(2)));
	  	  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(3)));
	  	  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(4)));
	  	  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(5)));
	  	}
	  	float pt_mean = (vmsa::pt_lowKS[energy][i_pt]+vmsa::pt_upKS[energy][i_pt])/2.0;

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
        }
      }
    }
  }
  
/*#if _PlotQA_
  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; ++i_theta)
  {
    c_diff->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_0_%s_Norm_%d",vmsa::pt_QAKS[energy],9,i_theta,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
    f_bw->SetParameter(0,Par[KEY_QA][0]);
    f_bw->SetParameter(1,Par[KEY_QA][1]);
    f_bw->SetParameter(2,Par[KEY_QA][2]);
    f_bw->SetLineColor(4);
    f_bw->SetLineStyle(2);
    f_bw->SetLineWidth(2);
    f_bw->Draw("l same");
  }

  c_diff->cd(9);
  for(int i_sigma = Sig_start; i_sigma < Sig_stop; ++i_sigma)
  {
    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
    {
      string KEY_counts_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_0_%s_Norm_%d_Sigma_%d_%s",vmsa::pt_QAKS[energy],9,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,i_sigma,vmsa::mInteMethod[i_method].c_str());
      if(i_sigma == Sig_start && i_method == vmsa::Method_start)
      {
	h_mCounts[KEY_counts_QA];
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
      }
      else
      {
	h_mCounts[KEY_counts_QA]->SetMarkerStyle(24);
	h_mCounts[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
	h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
	h_mCounts[KEY_counts_QA]->Draw("pE same");
      }
      TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
      f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
      f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
      f_rho->SetLineColor(i_sigma+10*i_method+1);
      f_rho->SetLineWidth(2);
      f_rho->SetLineStyle(2);
      f_rho->Draw("l same");
    }
  }
  if(!random3D) c_diff->SaveAs("../figures/c_diff_2.pdf");
  if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2.pdf");
#endif
*/
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

  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < 10/*vmsa::Cent_stop*/; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++) // systematic loop for nSigma
        {
          if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	  {
	    for(int i_sigma = Sig_start; i_sigma < Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
	        string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_0_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,0,1.1);
	      }
	    }
	  }
        }
      }
    }
  }
  if(!random3D) c_rho->SaveAs("../figures/c_rho.pdf");
  if(random3D)  c_rho->SaveAs("../figures/3DRandom/c_rho.pdf");

  string outputfile = Form("../output/AuAu%s/%s/RawKStarPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/RawKStarPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < 10/*vmsa::Cent_stop*/; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++) // systematic loop for nSigma
        {
          if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	  {
	    for(int i_sigma = Sig_start; i_sigma < Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
	        string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
	        g_mRho[KEY_rho]->Write();
	        for(int i_pt = vmsa::pt_rebin_firstKS[energy]; i_pt < vmsa::pt_rebin_lastKS[energy]; ++i_pt) // pt loop
	        {
	  	string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	  	h_mCounts[KEY_counts]->Write();
	        }
	      }
	    }
	  }
        }
      }
    }
  }
  File_OutPut->Close();
}
