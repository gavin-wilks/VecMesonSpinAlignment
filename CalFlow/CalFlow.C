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

void CalFlow(int energy = 4, int pid = 0, int year = 0)
{
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptDate(0);

  string inputfile = Form("../data/file_%s_Resolution.root",vmsa::mBeamEnergy[energy].c_str());
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

  string InPutFile_SE = "../data/FlowYields_Phi_SE_19GeV.root";
  
  string InPutFile_ME = "../data/FlowYields_Phi_ME_19GeV.root";
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
    {
      //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      //{
	for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
	{
	  string KEY_SE = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  string KEY_ME = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());

          if(i_cent == 0)
          {       
	    string KEY_SE10 = Form("pt_%d_Centrality_10_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",i_pt,i_theta,vmsa::mPID[pid].c_str());
	    h_mMass_SE[KEY_SE10] = (TH1F*)((TH1F*)File_SE->Get(KEY_SE.c_str()))->Clone(KEY_SE10.c_str()); 

	    string KEY_ME10 = Form("pt_%d_Centrality_10_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",i_pt,i_theta,vmsa::mPID[pid].c_str());
	    h_mMass_ME[KEY_ME10] = (TH1F*)((TH1F*)File_ME->Get(KEY_ME.c_str()))->Clone(KEY_ME10.c_str()); 
          }
          if(i_cent > 0 && i_cent < 10)
          {
	    string KEY_SE10 = Form("pt_%d_Centrality_10_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",i_pt,i_theta,vmsa::mPID[pid].c_str());
            h_mMass_SE[KEY_SE10]->Add((TH1F*)File_SE->Get(KEY_SE.c_str())->Clone());
	    string KEY_ME10 = Form("pt_%d_Centrality_10_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",i_pt,i_theta,vmsa::mPID[pid].c_str());
            h_mMass_ME[KEY_ME10]->Add((TH1F*)File_ME->Get(KEY_ME.c_str())->Clone());
            //h_mMass_SE[KEY_SE10]->Print();
            //h_mMass_ME[KEY_ME10]->Print();
          }
          if(i_cent == 10)
          {
            string KEY_SE = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][0]);
	    int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][0]);
	    float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);

	    string KEY_ME = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
            //if(isnan(Inte_ME)) Inte_ME = 0.0;
	    h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
            
            //std::cout << "Inte_SE = " << Inte_SE << ",    Inte_ME = " << Inte_ME << std::endl;

	    string KEY_SM = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	    h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
            continue;
          }
             
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
          
          
          //if(i_cent == 9) std::cout << "Inte_SE = " << Inte_SE << ",    Inte_ME = " << Inte_ME << std::endl;
          
	}
      //}
    }
  }

#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *cy6 = new TCanvas("cy6","cy6",10,10,800,600);
  cy6->cd();
  string KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",6,10,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",6,10,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",6,10,3,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy6->SaveAs("./figures/cy6_costheta3_pt6_cent9.pdf");

  TCanvas *cy7 = new TCanvas("cy7","cy7",10,10,800,600);
  cy7->cd();
  KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",7,10,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",7,10,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",7,10,3,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy7->SaveAs("./figures/cy7_costheta3_pt7_cent9.pdf");

  TCanvas *cy5 = new TCanvas("cy5","cy5",10,10,800,600);
  cy5->cd();
  KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",5,10,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",5,10,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",5,10,0,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy5->SaveAs("./figures/cy5_costheta3_pt5_cent9.pdf");
 


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
    string KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_SE",i_pt,10,3,vmsa::mPID[pid].c_str());
    h_mMass_SE[KEY_SE_QA]->DrawCopy();

    string KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_0_Sig_0_%s_ME",i_pt,10,3,vmsa::mPID[pid].c_str());
    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,10,3,vmsa::mPID[pid].c_str());
    h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
  {
    //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    //{
      for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
      {
	for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",pt_bin,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	  {
	    string KEY_SM = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    //	      cout << "KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	    if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
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
    //}
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
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,10,3,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetMarkerStyle(20);
    h_mMass[KEY_QA]->SetMarkerSize(0.4);
    h_mMass[KEY_QA]->SetLineColor(1);
    h_mMass[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  if(pid == 0) // Polynomial fit subtraction is only needed for phi meson
  {
    // Poly + Breit Wignar fit to phi integrated InvMass
    TH1FMap h_mMass_theta;
    vecFMap ParFit_theta;

    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
      {
	//for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
	//{
	  string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	  for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    if(i_theta == 0) h_mMass_theta[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone();
	    else h_mMass_theta[KEY_theta]->Add(h_mMass[KEY],1.0);
	  }
	  TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
	  for(int i_par = 0; i_par < 5; i_par++)
	  {
	    f_bw->ReleaseParameter(i_par);
	  }
	  f_bw->SetParameter(0,1.019);
	  f_bw->SetParLimits(0,1.014,1.024);
	  f_bw->SetParameter(1,0.0055);
	  f_bw->SetParameter(2,10000);
	  f_bw->SetParameter(3,-6000);
	  f_bw->SetParameter(4,0.5);
	  f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  ParFit_theta[KEY_theta].clear();
	  h_mMass_theta[KEY_theta]->Fit(f_bw,"NQR");
	  for(int n_par = 0; n_par < 5; n_par++)
	  {
	    ParFit_theta[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(n_par)));
	  }
	//}
      }
    }

#if _PlotQA_
    // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
    TCanvas *c_linsub = new TCanvas("c_linsub","c_linsub",10,10,800,600);
    TCanvas *c_mMass_theta = new TCanvas("c_mMass_theta","c_mMass_theta",10,10,1400,1400);
    c_mMass_theta->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,10,vmsa::mPID[pid].c_str()); 
      h_mMass_theta[KEY_theta_QA]->SetMarkerColor(1);
      h_mMass_theta[KEY_theta_QA]->SetMarkerStyle(24);
      h_mMass_theta[KEY_theta_QA]->SetMarkerSize(0.8);
      h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");

      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
      for(int i_par = 0; i_par < 5; i_par++)
      {
	f_bw->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par]);
      }
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->DrawCopy("l same");

      TF1 *f_poly = new TF1("f_poly",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
      f_poly->SetParameter(0,ParFit_theta[KEY_theta_QA][3]);
      f_poly->SetParameter(1,ParFit_theta[KEY_theta_QA][4]);
      f_poly->SetLineColor(4);
      f_poly->SetLineStyle(2);
      f_poly->SetLineWidth(4);
      f_poly->DrawCopy("l same");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);

      if(i_pt == 2)
      {
        c_linsub->cd();
        h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");
        f_bw->DrawCopy("l same");
        f_poly->DrawCopy("l same");

        string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
        plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
        PlotLine(0.98,1.05,0.0,0.0,1,2,2);
        c_linsub->SaveAs("./figures/c_linsub.pdf");
      }
    }
    c_mMass_theta->SaveAs("./figures/c_mMass_theta.pdf");
#endif

    // Poly+bw fits for phi differential InvMass
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
      {
	//for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
	//{
	  string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	  for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5); //Poly+bw fits
	    for(int i_par = 0; i_par < 5; i_par++)
	    {
	      f_bw->ReleaseParameter(i_par);
	    }
	    f_bw->FixParameter(0,ParFit_theta[KEY_theta][0]);
	    f_bw->FixParameter(1,ParFit_theta[KEY_theta][1]);
	    f_bw->SetParameter(2,ParFit_theta[KEY_theta][2]/10.0);
	    f_bw->SetParameter(3,ParFit_theta[KEY_theta][3]/10.0);
	    f_bw->SetParameter(4,ParFit_theta[KEY_theta][4]/10.0);
	    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	    h_mMass[KEY]->Fit(f_bw,"NQR");

	    TF1 *f_poly = new TF1("f_poly",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
	    f_poly->FixParameter(0,f_bw->GetParameter(3));
	    f_poly->FixParameter(1,f_bw->GetParameter(4));

	    h_mMass[KEY]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
	  }
	//}
      }
    }

#if _PlotQA_
    // QA plots for phi differential InvMass after linear background subtraction
    TCanvas *c_mMass_phi_diff = new TCanvas("c_mMass_phi_diff","c_mMass_phi_diff",10,10,2500,2500);
    c_mMass_phi_diff->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,10,0,vmsa::mPID[pid].c_str());
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("PE");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }
    c_mMass_phi_diff->SaveAs("./figures/c_mMass_phi_diff.pdf");
#endif
  }

  TH1FMap h_mMass_total; // cos(theta*) integrated InvMass after linear background subtraction for bw fits to extract yields 
  vecFMap ParBW;

  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
    {
      //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      //{
	string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  if(i_theta == 0) h_mMass_total[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone();
	  else h_mMass_total[KEY_theta]->Add(h_mMass[KEY],1.0);
	}
	TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
	f_bw->SetParameter(0,vmsa::InvMass[pid]);
	f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.001,vmsa::InvMass[pid]+0.001);
	f_bw->SetParameter(1,vmsa::Width[pid]);
	f_bw->SetParameter(2,1000);
	f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	h_mMass_total[KEY_theta]->Fit(f_bw,"MQNR");
	ParBW[KEY_theta].clear();
	ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(0)));
	ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(1)));
	ParBW[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(2)));
      }
    //}
  }

#if _PlotQA_
  // QA: bw fits to phi integrated InvMass
  TCanvas *c_mMass_bw = new TCanvas("c_mMass_bw","c_mMass_bw",10,10,1400,1400);
  c_mMass_bw->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,10,vmsa::mPID[pid].c_str());
    h_mMass_total[KEY_theta_QA]->SetMarkerColor(1);
    h_mMass_total[KEY_theta_QA]->SetMarkerStyle(24);
    h_mMass_total[KEY_theta_QA]->SetMarkerSize(0.8);
    h_mMass_total[KEY_theta_QA]->DrawCopy("PE");
    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
    f_bw->SetParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->SetParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
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
  TH1FMap h_mCounts, h_mv2;
  vecFMap ParSpin_Count, ParSpin_BW;
  for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
  {
    //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    //{
      string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      h_mv2[KEY_v2_Count] = new TH1F(KEY_v2_Count.c_str(),KEY_v2_Count.c_str(),100,-0.05,9.95);
      string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
      h_mv2[KEY_v2_BW] = new TH1F(KEY_v2_BW.c_str(),KEY_v2_BW.c_str(),100,-0.05,9.95);
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
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
	  h_mCounts[KEY_Count]->SetBinContent(h_mCounts[KEY_Count]->FindBin(bin_center),counts);
	  h_mCounts[KEY_Count]->SetBinError(h_mCounts[KEY_Count]->FindBin(bin_center),TMath::Sqrt(errors));

	  // breit wigner integrating
	  TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
	  f_bw->FixParameter(0,ParBW[KEY_theta][0]);
	  f_bw->FixParameter(1,ParBW[KEY_theta][1]);
	  f_bw->SetParameter(2,ParBW[KEY_theta][2]/10.0);
	  f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  h_mMass[KEY]->Fit(f_bw,"NMQRI");
	  float bin_width = h_mMass[KEY]->GetBinWidth(1);
	  float Inte_start = ParBW[KEY_theta][0]-vmsa::nSigVec*ParBW[KEY_theta][1]-0.5*bin_width;
	  float Inte_stop  = ParBW[KEY_theta][0]+vmsa::nSigVec*ParBW[KEY_theta][1]+0.5*bin_width;
	  float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
	  float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	  h_mCounts[KEY_BW]->SetBinContent(h_mCounts[KEY_BW]->FindBin(bin_center),counts_bw);
	  h_mCounts[KEY_BW]->SetBinError(h_mCounts[KEY_BW]->FindBin(bin_center),errors_bw); 
          
          if(i_cent < 9 && i_pt < 5)
          {
            //if(isnan(counts) || counts < 0.0) counts = 0.0;        
            //if(isnan(counts_bw) || counts_bw < 0.0) counts_bw = 0.0;        

            std::cout << "pt bin = " << i_pt << ",    i_cent = " << i_cent << std::endl;
            phiMesons[0][i_pt][i_cent] += counts;
            phiMesons[1][i_pt][i_cent] += counts_bw;
            std::cout << "IndividualBin: " << phiMesons[0][i_pt][i_cent] << std::endl;
            std::cout << "IndividualBin: " << phiMesons[1][i_pt][i_cent] << std::endl;

            phiMesons[0][i_pt][9] += counts;
            phiMesons[1][i_pt][9] += counts_bw;
            std::cout << "AllCent: " << phiMesons[0][i_pt][9] << std::endl;
            std::cout << "AllCent: " << phiMesons[1][i_pt][9] << std::endl;

            phiMesons[0][6][9] += counts;
            phiMesons[1][6][9] += counts_bw;
            std::cout << "All: " << phiMesons[0][6][9] << std::endl;
            std::cout << "All: " << phiMesons[1][6][9] << std::endl;
          } 
	}
	float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

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
      }
    //}
  }

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",2,10,vmsa::mPID[pid].c_str());
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
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",2,10,i_theta,vmsa::mPID[pid].c_str());
    // h_mMass[KEY_QA]->SetTitle(Title_PhiPsi[i_theta-2].c_str());
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

    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
    f_bw->FixParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->FixParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]/10.0);
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
  string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",2,10,vmsa::mPID[pid].c_str());
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

  string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",2,10,vmsa::mPID[pid].c_str());
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
  //string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  TCanvas *c_mMass_psi_2 = new TCanvas("c_mMass_psi_2","c_mMass_psi_2",10,10,1000,1400);
  c_mMass_psi_2->Divide(3,4);
  for(int i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi_2->cd(i_theta+1)->SetLeftMargin(0.15);
    c_mMass_psi_2->cd(i_theta+1)->SetBottomMargin(0.15);
    c_mMass_psi_2->cd(i_theta+1)->SetTicks(1,1);
    c_mMass_psi_2->cd(i_theta+1)->SetGrid(0,0);
    c_mMass_psi_2->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_%s_SM",2,10,i_theta,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetTitle(Form("#phi-Psi_{2} bin %d",i_theta));
    h_mMass[KEY_QA]->SetMarkerColor(1);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(0.8);
    h_mMass[KEY_QA]->DrawCopy("pE");

    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
    f_bw->FixParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->FixParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]/10.0);
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
  TF1 *f_bw_QA = new TF1("f_bw_QA",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
  f_bw_QA->FixParameter(0,ParBW[KEY_theta_QA][0]);
  f_bw_QA->FixParameter(1,ParBW[KEY_theta_QA][1]);
  f_bw_QA->FixParameter(2,ParBW[KEY_theta_QA][2]);
  f_bw_QA->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
  f_bw_QA->SetLineColor(2);
  f_bw_QA->Draw("l same");
  string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][2],vmsa::pt_up[energy][2]);
  plotTopLegend((char*)pT_range.c_str(),0.15,0.7,0.08,1,0.0,42,1);
  PlotLine(0.98,1.05,0.0,0.0,1,2,2);


  c_mMass_psi_2->cd(12);
  //string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
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

  //string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
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

  /*for(int i_cent = vmsa::Cent_start; i_cent < 9; i_cent++) // Centrality loop
  {
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++)
    {
      float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

      string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits

      float count_content = h_mv2[KEY_v2_Count]->GetBinContent(h_mv2[KEY_v2_Count]->FindBin(pt_mean)); // bin counting
      float count_error   = h_mv2[KEY_v2_Count]->GetBinError(h_mv2[KEY_v2_Count]->FindBin(pt_mean));

      float bw_content = h_mv2[KEY_v2_BW]->GetBinContent(h_mv2[KEY_v2_BW]->FindBin(pt_mean)); // breit wigner fits
      float bw_error   = h_mv2[KEY_v2_BW]->GetBinError(h_mv2[KEY_v2_BW]->FindBin(pt_mean));

      if(i_cent < 9)
      {
        phiMesons[0][i_pt][i_cent] += count_content;
        phiMesons[1][i_pt][i_cent] += bw_content;
        std::cout << "IndividualBin: " << phiMesons[0][i_pt][i_cent] << std::endl;
        std::cout << "IndividualBin: " << phiMesons[1][i_pt][i_cent] << std::endl;

        phiMesons[0][i_pt][9] += count_content;
        phiMesons[1][i_pt][9] += bw_content;
        std::cout << "AllCent: " << phiMesons[0][i_pt][i_cent] << std::endl;
        std::cout << "AllCent: " << phiMesons[1][i_pt][i_cent] << std::endl;

        phiMesons[0][6][9] += count_content;
        phiMesons[1][6][9] += bw_content;
        std::cout << "All: " << phiMesons[0][i_pt][i_cent] << std::endl;
        std::cout << "All: " << phiMesons[1][i_pt][i_cent] << std::endl;
      } 
    }
  }*/

  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < 5/*vmsa::pt_rebin_last[energy]*/; i_pt++)
  {
    resWeight[0][i_pt] = 0.0;
    resWeight[1][i_pt] = 0.0;
    for(int i_cent = vmsa::Cent_start; i_cent < 9; i_cent++) // Centrality loop
    {
      resWeight[0][i_pt] += (phiMesons[0][i_pt][i_cent] * (1.0/mTpcSubRes2Val[i_cent]));
      resWeight[1][i_pt] += (phiMesons[1][i_pt][i_cent] * (1.0/mTpcSubRes2Val[i_cent]));
      //resWeightI[0][i_pt][i_cent] = (phiMesons[0][i_pt][i_cent]/phiMesons[0][i_pt][9]) * (1.0/mTpcSubRes2Val[i_cent]);
    }
    resWeight[0][i_pt] /= phiMesons[0][i_pt][9];
    resWeight[1][i_pt] /= phiMesons[1][i_pt][9];
  }

  TGraMap g_mv2;
  for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
  {
    string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
    g_mv2[KEY_v2_Count] = new TGraphAsymmErrors();
    string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
    g_mv2[KEY_v2_BW] = new TGraphAsymmErrors();
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < 5/*vmsa::pt_rebin_last[energy]*/; i_pt++)
    {
      float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

      float count_content = h_mv2[KEY_v2_Count]->GetBinContent(h_mv2[KEY_v2_Count]->FindBin(pt_mean)); // bin counting
      float count_error   = h_mv2[KEY_v2_Count]->GetBinError(h_mv2[KEY_v2_Count]->FindBin(pt_mean));
      if(i_cent != 10)
      {
        g_mv2[KEY_v2_Count]->SetPoint(i_pt,pt_mean+0.05,count_content);
        g_mv2[KEY_v2_Count]->SetPointError(i_pt,0.0,0.0,count_error,count_error);
      }
      if(i_cent == 10)
      {
        g_mv2[KEY_v2_Count]->SetPoint(i_pt,pt_mean+0.05,count_content*resWeight[0][i_pt]);
        g_mv2[KEY_v2_Count]->SetPointError(i_pt,0.0,0.0,count_error*resWeight[0][i_pt],count_error*resWeight[0][i_pt]);
      }
      string Name_count = Form("v2_%s_Centrality_%d_Count",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
      g_mv2[KEY_v2_Count]->SetName(Name_count.c_str());
      g_mv2[KEY_v2_Count]->SetMarkerStyle(24);
      g_mv2[KEY_v2_Count]->SetMarkerColor(4);

      float bw_content = h_mv2[KEY_v2_BW]->GetBinContent(h_mv2[KEY_v2_BW]->FindBin(pt_mean)); // breit wigner fits
      float bw_error   = h_mv2[KEY_v2_BW]->GetBinError(h_mv2[KEY_v2_BW]->FindBin(pt_mean));
      if(i_cent != 10)
      {
        g_mv2[KEY_v2_BW]->SetPoint(i_pt,pt_mean,bw_content);
        g_mv2[KEY_v2_BW]->SetPointError(i_pt,0.0,0.0,bw_error,bw_error);
      }
      if(i_cent == 10)
      {
        g_mv2[KEY_v2_BW]->SetPoint(i_pt,pt_mean,bw_content*resWeight[1][i_pt]);
        g_mv2[KEY_v2_BW]->SetPointError(i_pt,0.0,0.0,bw_error*resWeight[1][i_pt],bw_error*resWeight[1][i_pt]);
      }
      string Name_bw = Form("v2_%s_Centrality_%d_BW",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
      g_mv2[KEY_v2_BW]->SetName(Name_bw.c_str());
      g_mv2[KEY_v2_BW]->SetMarkerStyle(24);
      g_mv2[KEY_v2_BW]->SetMarkerColor(2);
    }
  }

  //std::cout << "Weighted Mean (Count)    = " << weightedMeanC/weightedErrorC << " +/- " << TMath::Sqrt(1.0/weightedErrorC) << std::endl;
  //std::cout << "Weighted Mean (Integral) = " << weightedMeanI/weightedErrorI << " +/- " << TMath::Sqrt(1.0/weightedErrorI) << std::endl;
  
  TFile *besi = TFile::Open("../data/BESI/HEPData-ins1210464-v1-root.root");
  TDirectory *dir = (TDirectory*) besi->Get("Table 200;1");
  dir->cd();
  TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1;1");
  // QA Plots for rho00 vs. pt of bin counting and breit wigner integrating
  TCanvas *c_v2_pT = new TCanvas("c_v2_pT","c_v2_pT",10,10,800,800);
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

  string KEY_v2_Count_QA = Form("v2_Count_Centrality_%d_%s",10,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2[KEY_v2_Count_QA] ,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.28,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.275,0.03,1,0.0,42,0);

  string KEY_v2_BW_QA = Form("v2_BW_Centrality_%d_%s",10,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mv2[KEY_v2_BW_QA] ,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.255,0.03,1,0.0,42,0);

  string leg_besi = "BES-I";
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19,25,1,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,25,1,1.3);
  plotTopLegend((char*)leg_besi.c_str(),0.6,0.235,0.03,1,0.0,42,0);

  //PlotLine(3.0,3.5,0.45,0.45,1,2,2);
  //string leg_line = "#rho_{00} = 1/3";
  //plotTopLegend((char*)leg_line.c_str(),3.6,0.447,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.275,0.04,1,0.0,42,0);
  string leg_centrality = "0%-80%";
  plotTopLegend((char*)leg_centrality.c_str(),3.1,0.245,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/v2_pT_%s.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_v2_pT->SaveAs(figure_name.c_str());

  
  string OutPutFile = "../output/RawPhiv2/RawPhiV2PtSys.root"; 
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_play->Write();
  for(int i_cent = vmsa::Cent_start; i_cent < 11; i_cent++) // Centrality loop
  {
    //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    //{
      string KEY_v2_Count = Form("v2_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      g_mv2[KEY_v2_Count]->Write();
      string KEY_v2_BW = Form("v2_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
      g_mv2[KEY_v2_BW]->Write();
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
      {
	string KEY_Count = Form("v2_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_Count]->Write();
	string KEY_BW = Form("v2_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_BW]->Write();
      }
    //}
  }
  File_OutPut->Close();
}
