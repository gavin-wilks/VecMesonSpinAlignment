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
#include "TFitResultPtr.h" 

#ifndef _PlotQA_
#define _PlotQA_  1
#endif


void calFlowSysKStar_Poly3(int energy = 4, int pid = 2, int year = 0, int ptQA = 5, int centQA = 6, int dcaQA = 0, int nsigQA = 0, int nhitQA = 0)
{
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0000001);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);
  string inputfile = Form("../output/AuAu%s/Flow/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
	  {
            if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
	    for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	    {
	      string KEY_InteTheta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              //int thetabin = 0;
	      for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++) // cos(theta*) loop
	      {
	        string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        /*if(i_theta%2 == 0)*/ h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str())->Clone();
	        //if(i_theta%2 == 1) 
                //{
                //  TH1F* hMassTemp = (TH1F*)File_InPut->Get(KEY.c_str())->Clone(); 
                //  h_mMass[KEY]->Add(hMassTemp,1.0);
                  //thetabin++;
                //}
                cout << "Loading --> " << KEY << endl;
	        //h_mMass[KEY]->Rebin(2);

	        if(i_theta == vmsa::PhiPsi_start) h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)h_mMass[KEY]->Clone(KEY_InteTheta.c_str());
	        else h_mMass_InteTheta[KEY_InteTheta]->Add(h_mMass[KEY],1.0);
	      }
	      TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	      f_bw->SetParameter(0,vmsa::InvMass[pid]);
	      f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.03,vmsa::InvMass[pid]+0.03);
	      f_bw->SetParameter(1,vmsa::Width[pid]);
	      f_bw->SetParameter(2,1.0);
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
	      Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(6)));
              delete f_bw;
	    }
          }
	}
      }
    }
  }



/*#if _PlotQA_
  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,1200);
  c_diff->Divide(3,4);
  for(int i_theta = 0; i_theta < 12; i_theta++)
  {
    c_diff->cd(i_theta+1);
    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
    c_diff->cd(i_theta+1)->SetTicks(1,1);
    c_diff->cd(i_theta+1)->SetGrid(0,0);
    if(i_theta < vmsa::PhiPsi_stop)
    {
      string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",ptQA,centQA,i_theta,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass[KEY_QA]->SetTitle("");
      h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
      h_mMass[KEY_QA]->GetXaxis()->SetRangeUser(0.75,1.08);

      h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
      h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
      h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
      h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
      h_mMass[KEY_QA]->GetYaxis()->CenterTitle();

      //h_mMass[KEY_QA]->GetXaxis()->SetRangeUser(0.75,1.1);

      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
      h_mMass[KEY_QA]->SetMarkerSize(1.2);
      h_mMass[KEY_QA]->Draw("pE");
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
    if(i_theta == vmsa::PhiPsi_stop)
    {
      string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",ptQA,centQA,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetRangeUser(0.75,1.08);

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
      TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][0]);
      f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][1]);
      f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][2]);
      f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA][3]);
      f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA][4]);
      f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA][5]);
      f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_QA][6]);
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(2);
      f_bw->SetLineWidth(2);
      f_bw->Draw("l same");

      TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
      f_bg->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][3]);
      f_bg->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][4]);
      f_bg->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][5]);
      f_bg->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA][6]);
      f_bg->SetLineColor(4);
      f_bg->SetLineStyle(2);
      f_bg->SetLineWidth(2);
      f_bg->Draw("l same");
    }
  }
#endif
*/
  // extract counts vs. pT with diffenretial integration ranges and methods
  vecFMap Par;
  TH1FMap h_mCounts;
  vecFMap Par_v2Fit;
  TGraMap g_mV2;
  for(int i_cent = vmsa::Cent_start; i_cent < 13; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
	{
          if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	  {
	    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
	        string KEY_v2 = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        g_mV2[KEY_v2] = new TGraphAsymmErrors();
                
                float weight;
                string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
                string KEY_WeightBW   = Form("ResWeight_BW_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm); 

	        for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; ++i_pt) // pt loop
	        {
	  	  string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	  	  h_mCounts[KEY_counts] = new TH1F(KEY_counts.c_str(),KEY_counts.c_str(),10,0.0,TMath::Pi()/2.0);

	  	  string KEY_InteTheta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	  	  for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++) // cos(theta*) loop
	  	  {
	  	    string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
                    cout << KEY << endl;
	  	    TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	  	    f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta][0]);
	  	    f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta][1]);
	  	    f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta][2]/10.0);
	  	    f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta][3]/10.0);
	  	    f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta][4]/10.0);
	  	    f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta][5]);
                    f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta][6]);
	  	    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  	    TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"MQNRS");
	  	    float bin_center = (vmsa::PhiPsi_low[i_theta] + vmsa::PhiPsi_up[i_theta])/2.0;
                    
                    TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
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
	  	    float Inte_start = Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSysKS[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
                    float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSysKS[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
                    float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
                    float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width;

	  	    if(i_method == 0)
	  	    {
	  	      int bin_start = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSysKS[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
	  	      int bin_stop  = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSysKS[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
	  	      float counts = 0.0;
	  	      float errors = 0.0;
	  	      for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	  	      {
	  	        counts += h_mMass[KEY]->GetBinContent(i_bin);
	  	        errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	  	      }
	  	      h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts);
	  	      h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors));
                      std::vector<float> *v_weight; 
                      File_InPut->GetObject(KEY_WeightGaus.c_str(),v_weight);
                      weight = v_weight->at(0);
                      cout << "Loading --> " << KEY_WeightGaus << endl;
                      cout << "weight = " << weight << endl;
	  	    }
	  	    if(i_method == 1)
	  	    {
	  	      float bin_width = h_mMass[KEY]->GetBinWidth(1);
	  	      float Inte_start = Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSysKS[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
	  	      float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSysKS[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
	  	      float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
	  	      float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	  	      h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts_bw);
	  	      h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),errors_bw);
                      std::vector<float> *v_weight; 
                      File_InPut->GetObject(KEY_WeightBW.c_str(),v_weight);
                      weight = v_weight->at(0);
                      cout << "Loading --> " << KEY_WeightBW << endl;
                      cout << "weight = " << weight << endl;
	  	    }
	  	    Par[KEY].clear();
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(0)));
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(1)));
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(2)));
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(3)));
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(4)));
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(5)));
	  	    Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(6)));
                    delete f_bw;
                    delete f_bg;
	  	  }
	  	  float pt_mean = (vmsa::pt_lowKSv2[energy][i_pt]+vmsa::pt_upKSv2[energy][i_pt])/2.0;

	  	  TF1 *f_v2 = new TF1("f_v2",flow,0.0,TMath::Pi()/2.0,2);
	  	  f_v2->SetParameter(0,0.05);
	  	  f_v2->SetParameter(1,h_mCounts[KEY_counts]->GetMaximum());
	  	  h_mCounts[KEY_counts]->Fit(f_v2,"NMRI");
	  	  Par_v2Fit[KEY_counts].clear();
	  	  Par_v2Fit[KEY_counts].push_back(static_cast<float>(f_v2->GetParameter(0)));
	  	  Par_v2Fit[KEY_counts].push_back(static_cast<float>(f_v2->GetParameter(1)));

	  	  g_mV2[KEY_v2]->SetPoint(i_pt,pt_mean,f_v2->GetParameter(0)*weight);
	  	  g_mV2[KEY_v2]->SetPointError(i_pt,0.0,0.0,f_v2->GetParError(0)*weight,f_v2->GetParError(0)*weight); 
	        }
	      }
	    }
	  }
        }
      }
    }
  }
  
#if _PlotQA_

  for(int i_cent = vmsa::Cent_start; i_cent < 13; ++i_cent) // Centrality loop
  {
    string outputname = Form("./figures/Phi/19GeV/allThetaYieldsPoly3_cent%d.pdf",i_cent);
    string output_start = Form("%s[",outputname.c_str());
    
    TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,1200);
    c_diff->Divide(3,4);
    c_diff->Print(output_start.c_str());
    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; ++i_pt) // pt loop
    {

      for(int i_theta = 0; i_theta < 12; i_theta++)
      {
        c_diff->cd(i_theta+1);
        c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
        c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
        c_diff->cd(i_theta+1)->SetTicks(1,1);
        c_diff->cd(i_theta+1)->SetGrid(0,0);
        if(i_theta < vmsa::PhiPsi_stop)
        {
          string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
          h_mMass[KEY_QA]->SetTitle("");
          h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
          h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
          h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
          h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
          h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
          h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
          h_mMass[KEY_QA]->GetXaxis()->SetRangeUser(0.75,1.08);

          h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
          h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
          h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
          h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
          h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
          h_mMass[KEY_QA]->GetYaxis()->CenterTitle();

          //h_mMass[KEY_QA]->GetXaxis()->SetRangeUser(0.75,1.1);

          h_mMass[KEY_QA]->SetMarkerStyle(24);
          h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
          h_mMass[KEY_QA]->SetMarkerSize(1.2);
          h_mMass[KEY_QA]->Draw("pE");
          PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
        }
        if(i_theta == vmsa::PhiPsi_stop)
        {
          string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
          h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
          h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetRangeUser(0.75,1.08);

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
          TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
          f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][0]);
          f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][1]);
          f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][2]);
          f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA][3]);
          f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA][4]);
          f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA][5]);
          f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_QA][6]);
          f_bw->SetLineColor(2);
          f_bw->SetLineStyle(2);
          f_bw->SetLineWidth(2);
          f_bw->Draw("l same");

          TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
          f_bg->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][3]);
          f_bg->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][4]);
          f_bg->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][5]);
          f_bg->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA][6]);
          f_bg->SetLineColor(4);
          f_bg->SetLineStyle(2);
          f_bg->SetLineWidth(2);
          f_bg->Draw("l same");
        }
      }
      for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++)
      {
        c_diff->cd(i_theta+1);
        string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        string KEY_InteTheta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
        f_bw->SetParameter(0,Par[KEY_QA][0]);
        f_bw->SetParameter(1,Par[KEY_QA][1]);
        f_bw->SetParameter(2,Par[KEY_QA][2]);
        f_bw->SetParameter(3,Par[KEY_QA][3]);
        f_bw->SetParameter(4,Par[KEY_QA][4]);
        f_bw->SetParameter(5,Par[KEY_QA][5]);
        f_bw->SetParameter(6,Par[KEY_QA][6]);
        f_bw->SetLineColor(4);
        f_bw->SetLineStyle(2);
        f_bw->SetLineWidth(2);
        f_bw->Draw("l same");

        TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
        f_bg->SetParameter(0,Par[KEY_QA][3]);
        f_bg->SetParameter(1,Par[KEY_QA][4]);
        f_bg->SetParameter(2,Par[KEY_QA][5]);
        f_bg->SetParameter(3,Par[KEY_QA][6]);
        f_bg->SetLineColor(2);
        f_bg->SetLineStyle(2);
        f_bg->SetLineWidth(2);
        f_bg->Draw("l same");
        
        float x1_0 = Par_InteTheta[KEY_InteTheta][0] - 1.5*Par_InteTheta[KEY_InteTheta][1];
        float x2_0 = Par_InteTheta[KEY_InteTheta][0] + 1.5*Par_InteTheta[KEY_InteTheta][1];
        float x1_1 = Par_InteTheta[KEY_InteTheta][0] - 2.0*Par_InteTheta[KEY_InteTheta][1];
        float x2_1 = Par_InteTheta[KEY_InteTheta][0] + 2.0*Par_InteTheta[KEY_InteTheta][1];
        float x1_2 = Par_InteTheta[KEY_InteTheta][0] - 2.5*Par_InteTheta[KEY_InteTheta][1];
        float x2_2 = Par_InteTheta[KEY_InteTheta][0] + 2.5*Par_InteTheta[KEY_InteTheta][1];
        float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(Par_InteTheta[KEY_InteTheta][0]));
        //PlotLine(x1_0,x1_0,0,y,4,2,2);
        //PlotLine(x2_0,x2_0,0,y,4,2,2);
        PlotLine(x1_1,x1_1,0,y,1,2,2);
        PlotLine(x2_1,x2_1,0,y,1,2,2);
        //PlotLine(x1_2,x1_2,0,y,4,2,2);
        //PlotLine(x2_2,x2_2,0,y,4,2,2);
      }

      c_diff->cd(12);
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          string KEY_counts_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,i_sigma,vmsa::mInteMethod[i_method].c_str());
          if(i_sigma == vmsa::Sig_start && i_method == 1/*vmsa::Method_start*/)
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
          else if(i_method == 1)
          {
            h_mCounts[KEY_counts_QA]->SetMarkerStyle(24);
            h_mCounts[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
            h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
            h_mCounts[KEY_counts_QA]->Draw("pE same");
          }
          TF1 *f_v2 = new TF1("f_v2",flow,0.0,TMath::Pi()/2.0,2);
          f_v2->SetParameter(0,Par_v2Fit[KEY_counts_QA][0]);
          f_v2->SetParameter(1,Par_v2Fit[KEY_counts_QA][1]);
          f_v2->SetLineColor(i_sigma+10*i_method+1);
          f_v2->SetLineWidth(2);
          f_v2->SetLineStyle(2);
          f_v2->Draw("l same");
        }
      }
      c_diff->Update();
      c_diff->Print(outputname.c_str());
    }
    string output_stop = Form("%s]",outputname.c_str());
    c_diff->Print(output_stop.c_str()); // close pdf file
  }

#endif

  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd();
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
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

  h_frame->GetYaxis()->SetRangeUser(-0.1,0.3);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("v_{2}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  //PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  //for(int i_cent = vmsa::Cent_start; i_cent < 13; ++i_cent) // Centrality loop
  //{
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
	{
          if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	  {
	    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
	        string KEY_v2 = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",centQA,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mV2[KEY_v2],24,i_sigma+10*i_method+1,0,1.1);
	      }
	    }
	  }
        }
      }
    }
  //}
  c_v2->SaveAs("../figures/c_v2.pdf");

  string outputfile = Form("../output/AuAu%s/Flow/%s/RawKStarPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawV2PtSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_cent = vmsa::Cent_start; i_cent < 13; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
	{
          if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	  {
	    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
	        string KEY_v2 = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        g_mV2[KEY_v2]->SetName(KEY_v2.c_str());
	        g_mV2[KEY_v2]->Write();
                cout << "Writing --> " << KEY_v2 << endl;
	        for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; ++i_pt) // pt loop
	        {
	  	  string KEY_counts = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	  	  h_mCounts[KEY_counts]->Write();
                  cout << "Writing --> " << KEY_counts << endl;
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
