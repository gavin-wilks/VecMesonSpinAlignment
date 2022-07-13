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
#include "TFitResultPtr.h" 
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calSpinAlignmentKStarPoly3(int energy = 4, int pid = 2, int year = 0, double nSigVec = 2.0, int ptQA = 1)
{
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptDate(0);

  string InPutFile_SE = "../data/Yields_KStar_SE_19GeV_20220413.root";  
  string InPutFile_ME = "../data/Yields_KStar_ME_19GeV_20220429.root";

  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY_SE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  h_mMass_SE[KEY_SE] = (TH1F*)File_SE->Get(KEY_SE.c_str())->Clone(); 
          //h_mMass_SE[KEY_SE]->Rebin(4);
	  int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][0]);
	  int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][0]);
	  float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);

	  string KEY_ME = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  h_mMass_ME[KEY_ME] = (TH1F*)File_ME->Get(KEY_ME.c_str())->Clone(); 
          //h_mMass_ME[KEY_ME]->Rebin(4);
	  float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	  h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
	  string KEY_SM = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	  h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
	}
    }
  }

#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *cy6 = new TCanvas("cy6","cy6",10,10,800,600);
  cy6->cd();
  string KEY_SE_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",6,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",6,9,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",6,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy6->SaveAs("./figures/KStar/cy6_costheta3_pt6_cent9.pdf");

  TCanvas *cy7 = new TCanvas("cy7","cy7",10,10,800,600);
  cy7->cd();
  KEY_SE_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",7,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  KEY_ME_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",7,9,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  KEY_SM_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",7,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy7->SaveAs("./figures/KStar/cy7_costheta3_pt7_cent9.pdf");

  TCanvas *cy5 = new TCanvas("cy5","cy5",10,10,800,600);
  cy5->cd();
  KEY_SE_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",5,9,3,vmsa::mPID[pid].c_str());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  KEY_ME_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",5,9,3,vmsa::mPID[pid].c_str());
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  KEY_SM_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",5,9,vmsa::CTS_start,vmsa::mPID[pid].c_str());
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
  cy5->SaveAs("./figures/KStar/cy5_costheta3_pt5_cent9.pdf");
 


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
    string KEY_SE_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",i_pt,9,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    h_mMass_SE[KEY_SE_QA]->DrawCopy();

    string KEY_ME_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_ME",i_pt,9,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,9,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME (has residual background)
  TH1FMap h_mMassSBR; // rebinned InvMass distribution, SE (has residual background)
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    //{
      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
      {
	for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",pt_bin,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  string KEY_SBR = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",pt_bin,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	  {
	    string KEY_SM = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    string KEY_SE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
	    {
	      h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone();
              h_mMassSBR[KEY_SBR] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
            }
	    else
	    {
	      h_mMass[KEY]->Add(h_mMass_SM[KEY_SM],1.0);
	      h_mMassSBR[KEY_SBR]->Add(h_mMass_SE[KEY_SE],1.0);
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
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,9,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetMarkerStyle(20);
    h_mMass[KEY_QA]->SetMarkerSize(0.4);
    h_mMass[KEY_QA]->SetLineColor(1);
    h_mMass[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  //if(pid == 0 || pid == 2) // Polynomial fit subtraction is only needed for phi meson
  //{
    // Poly + Breit Wignar fit to phi integrated InvMass
    TH1FMap h_mMass_theta;
    vecFMap ParFit_theta;

    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
      {
	//for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
	//{
	  string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
            h_mMass[KEY]->Rebin(4);
	    if(i_theta == vmsa::CTS_start) h_mMass_theta[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone();
	    else h_mMass_theta[KEY_theta]->Add(h_mMass[KEY],1.0);
	  }
	  TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	  for(int i_par = 0; i_par < 7; i_par++)
	  {
	    f_bw->ReleaseParameter(i_par);
	  }
	  f_bw->SetParameter(0,vmsa::InvMass[pid]); //mass
	  f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.01,vmsa::InvMass[pid]+0.01); //mass limits
	  f_bw->SetParameter(1,vmsa::Width[pid]); //decay width
	  f_bw->SetParameter(2,10000); //normalization
	  //f_bw->SetParameter(3,0.0); //linear y-intercept
	  //f_bw->SetParameter(4,0.0); //linear slope
	  //f_bw->SetParameter(5,0.0); //linear slope
	  f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  ParFit_theta[KEY_theta].clear();
	  h_mMass_theta[KEY_theta]->Fit(f_bw,"NQR");
	  for(int n_par = 0; n_par < 7; n_par++)
	  {
	    ParFit_theta[KEY_theta].push_back(static_cast<float>(f_bw->GetParameter(n_par)));
            if(i_cent == 9) cout << "n_par = " << n_par << ", par = " << f_bw->GetParameter(n_par) << endl;
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

      TF1 *f_poly = new TF1("f_poly",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
      f_poly->SetParameter(0,ParFit_theta[KEY_theta_QA][3]);
      f_poly->SetParameter(1,ParFit_theta[KEY_theta_QA][4]);
      f_poly->SetParameter(2,ParFit_theta[KEY_theta_QA][5]);
      f_poly->SetParameter(3,ParFit_theta[KEY_theta_QA][6]);
      f_poly->SetLineColor(4);
      f_poly->SetLineStyle(2);
      f_poly->SetLineWidth(4);
      f_poly->DrawCopy("l same");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.6,1.2,0.0,0.0,1,2,2);

      if(i_pt == ptQA)
      {
        c_linsub->cd();
        h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");
        f_bw->DrawCopy("l same");
        f_poly->DrawCopy("l same");

        string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
        plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
        PlotLine(0.6,1.2,0.0,0.0,1,2,2);
        c_linsub->SaveAs("./figures/KStar/c_linsub.pdf");
      }
    }
    c_mMass_theta->SaveAs("./figures/KStar/c_mMass_theta.pdf");
#endif

    TH1FMap h_mMassBG;
    // Poly+bw fits for phi differential InvMass
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
      {
	//for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
	//{
	  string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	    //TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6); //Poly+bw fits
	    //for(int i_par = 0; i_par < 6; i_par++)
	    //{
	    //  f_bw->ReleaseParameter(i_par);
	    //}
	    //f_bw->FixParameter(0,ParFit_theta[KEY_theta][0]);
	    //f_bw->FixParameter(1,ParFit_theta[KEY_theta][1]);
	    //f_bw->SetParameter(2,ParFit_theta[KEY_theta][2]/7.0);
	    //f_bw->SetParameter(3,ParFit_theta[KEY_theta][3]/7.0);
	    //f_bw->SetParameter(4,ParFit_theta[KEY_theta][4]/7.0);
	    //f_bw->SetParameter(5,ParFit_theta[KEY_theta][5]/7.0);
	    //f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	    //TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"NQRS");

	    //TF1 *f_poly = new TF1("f_poly",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
	    //f_poly->SetParameter(0,f_bw->GetParameter(3));
	    //f_poly->SetParameter(1,f_bw->GetParameter(4));
	    //f_poly->SetParameter(2,f_bw->GetParameter(5));
	    //f_poly->SetParError(0,f_bw->GetParError(3));
	    //f_poly->SetParError(1,f_bw->GetParError(4));
	    //f_poly->SetParError(2,f_bw->GetParError(5));
           
            //double params[3] = {result->GetParams()[3],result->GetParams()[4],result->GetParams()[5]};
            //TMatrixDSym covArr(3);
            //covArr(0,0) = result->GetCovarianceMatrix()(3,3);
            //covArr(0,1) = result->GetCovarianceMatrix()(3,4);
            //covArr(0,2) = result->GetCovarianceMatrix()(3,5);
            //covArr(1,0) = result->GetCovarianceMatrix()(4,3);
            //covArr(1,1) = result->GetCovarianceMatrix()(4,4);
            //covArr(1,2) = result->GetCovarianceMatrix()(4,5);
            //covArr(2,0) = result->GetCovarianceMatrix()(5,3);
            //covArr(2,1) = result->GetCovarianceMatrix()(5,4);
            //covArr(2,2) = result->GetCovarianceMatrix()(5,5);

            //if(i_cent == 9)
            //{
            //  //cout << "fit parameters directly from fit" << endl;
            //  //cout << "intercept: " << f_poly->GetParameter(0) << endl;
            //  //cout << "error in intercept: " << f_poly->GetParError(0) << endl;
            //  //cout << "slope:     " << f_poly->GetParameter(1) << endl;
            //  //cout << "error in slope:     " << f_poly->GetParError(1) << endl;
            //  //cout << "fit parameters from fit result (COV)" << endl;
            //  //cout << "intercept: " << result->GetParams()[3] << endl;
            //  //cout << "error in intercept: " << covArr[0] << endl;
            //  //cout << "slope:     " << result->GetParams()[4] << endl;
            //  //cout << "error in slope:     " << covArr[1] << endl;
            //  //for(int i=0; i<30; i++)  cout << "getmatrixarray: i = " <<i << "  value = "<< result->GetCovarianceMatrix().GetMatrixArray()[i]<< endl;
            //}
      
          
            //TH1F *hTotal = (TH1F*)h_mMass[KEY]->Clone();
            h_mMassBG[KEY] = (TH1F*)h_mMass[KEY]->Clone();
	    //h_mMass[KEY]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
            //for(int bin = 1; bin <= h_mMass[KEY]->GetNbinsX(); bin++)
            //{
            //  float binWidth = h_mMass[KEY]->GetBinWidth(bin);  
            //  float binLowEdge = hTotal->GetBinLowEdge(bin);
            //  float funcError = f_poly->IntegralError(binLowEdge,binLowEdge+binWidth,params,covArr.GetMatrixArray())/binWidth;
            //  float totError = hTotal->GetBinError(bin);
            //  if(i_cent == 9)
            //  { 
            //  //  cout << "binWidth  = " << binWidth << endl;
            //  //  cout << "binLow    = " << binLowEdge << endl;
            //    //cout << "funcError = " << funcError << endl;
            //    //cout << "totError  = " << totError << endl;
            //  }
            //  h_mMass[KEY]->SetBinError(bin,TMath::Sqrt(totError*totError+funcError*funcError)); 
            //} 
            //delete hTotal;
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
      string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,9,vmsa::CTS_start,vmsa::mPID[pid].c_str());
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("PE");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.6,1.2,0.0,0.0,1,2,2);
    }
    c_mMass_phi_diff->SaveAs("./figures/KStar/c_mMass_phi_diff.pdf");
#endif
  

  TH1FMap h_mMass_total; // cos(theta*) integrated InvMass after linear background subtraction for bw fits to extract yields 
  TH1FMap h_mMass_BGtotal; // cos(theta*) integrated InvMass after linear background subtraction for bw fits to extract yields 
  vecFMap ParBW;

  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
	string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  if(i_theta == vmsa::CTS_start){ h_mMass_total[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone(); h_mMass_BGtotal[KEY_theta] = (TH1F*)h_mMassBG[KEY]->Clone();}
	  else {h_mMass_total[KEY_theta]->Add(h_mMass[KEY],1.0); h_mMass_BGtotal[KEY_theta]->Add(h_mMassBG[KEY],1.0);}

	}
	TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	f_bw->SetParameter(0,vmsa::InvMass[pid]);
	f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.01,vmsa::InvMass[pid]+0.01);
	f_bw->SetParameter(1,vmsa::Width[pid]);
	f_bw->SetParameter(2,1000);
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
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_mMass_bw->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,9,vmsa::mPID[pid].c_str());
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
    PlotLine(0.6,1.2,0.0,0.0,1,2,2);
  }
#endif

  // calculate counts and errors for cos(theta*) bin with bin counting and integrating
  TH1FMap h_mCounts, h_mRho00;
  vecFMap ParSpin_Count, ParSpin_BW;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
      string KEY_Rho00_Count = Form("Rho00_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      h_mRho00[KEY_Rho00_Count] = new TH1F(KEY_Rho00_Count.c_str(),KEY_Rho00_Count.c_str(),100,-0.05,9.95);
      string KEY_Rho00_BW = Form("Rho00_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
      h_mRho00[KEY_Rho00_BW] = new TH1F(KEY_Rho00_BW.c_str(),KEY_Rho00_BW.c_str(),100,-0.05,9.95);
      double aveSBRbw = 0;
      double aveSBRcount = 0;
      int totalSBR = 0;
      double lowSigBW = 100000;
      double lowSigCount = 100000;
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
      {
	string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_Count] = new TH1F(KEY_Count.c_str(),KEY_Count.c_str(),7,0.0,1.0);
	string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_BW] = new TH1F(KEY_BW.c_str(),KEY_BW.c_str(),7,0.0,1.0);
	string KEY_theta = Form("pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	float signal_counts = 0.0;
	float signal_counts_bw = 0.0;
        float SBRcounts = 0.0;
        totalSBR++;
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  // bin counting
	  string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  string KEY_SE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_0_Sig_0_NHit_0_%s_SE",i_pt,i_cent,i_theta,vmsa::mPID[pid].c_str());
	  float counts = 0.0;
	  float errors = 0.0;
          //float SBRcounts = 0.0;
	  float bin_center = 1/14.0+i_theta/7.0;
	  int bin_start = h_mMass[KEY]->FindBin(ParBW[KEY_theta][0]-nSigVec*ParBW[KEY_theta][1]);
	  int bin_stop  = h_mMass[KEY]->FindBin(ParBW[KEY_theta][0]+nSigVec*ParBW[KEY_theta][1]);
      
          

	  for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	  {
	    counts += h_mMass[KEY]->GetBinContent(i_bin);
	    if(i_cent == 9)
            {
              signal_counts += h_mMass[KEY]->GetBinContent(i_bin);
	      SBRcounts += h_mMassSBR[KEY_SE]->GetBinContent(i_bin);
              cout <<  h_mMass[KEY]->GetBinContent(i_bin) << endl;
            }
	    errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	  }
	  h_mCounts[KEY_Count]->SetBinContent(h_mCounts[KEY_Count]->FindBin(bin_center),counts);
	  h_mCounts[KEY_Count]->SetBinError(h_mCounts[KEY_Count]->FindBin(bin_center),TMath::Sqrt(errors));

	  // breit wigner integrating
	  TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	  f_bw->FixParameter(0,ParBW[KEY_theta][0]);
	  f_bw->FixParameter(1,ParBW[KEY_theta][1]);
	  f_bw->SetParameter(2,ParBW[KEY_theta][2]/7.0);
	  f_bw->SetParameter(3,ParBW[KEY_theta][3]/7.0);
	  f_bw->SetParameter(4,ParBW[KEY_theta][4]/7.0);
	  f_bw->SetParameter(5,ParBW[KEY_theta][5]/7.0);
	  f_bw->SetParameter(6,ParBW[KEY_theta][6]/7.0);
	  f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"NMQRIS");
	  float bin_width = h_mMass[KEY]->GetBinWidth(1);
	  float Inte_start = ParBW[KEY_theta][0]-nSigVec*ParBW[KEY_theta][1]-0.5*bin_width;
	  float Inte_stop  = ParBW[KEY_theta][0]+nSigVec*ParBW[KEY_theta][1]+0.5*bin_width;
	  float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
          if(i_cent == 9) signal_counts_bw+=counts_bw;
	  float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;

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

          float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
          float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width;
	  h_mCounts[KEY_BW]->SetBinContent(h_mCounts[KEY_BW]->FindBin(bin_center),counts_bw-counts_bg);
	  h_mCounts[KEY_BW]->SetBinError(h_mCounts[KEY_BW]->FindBin(bin_center),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));
	  h_mCounts[KEY_Count]->SetBinContent(h_mCounts[KEY_Count]->FindBin(bin_center),counts-counts_bg);
	  h_mCounts[KEY_Count]->SetBinError(h_mCounts[KEY_Count]->FindBin(bin_center),TMath::Sqrt(errors+errors_bg*errors_bg));
	}
	float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;
   
        if (i_cent == 9)
        {
          cout << "COUNTS = " << signal_counts << " | SBR-R COUNTS = " << SBRcounts << endl;
          cout << "PT = " << i_pt << " | COUNTING SIGNIFICANCE = " << signal_counts/TMath::Sqrt(SBRcounts) << endl;
          cout << "PT = " << i_pt << " | INTEGRAL SIGNIFICANCE = " << signal_counts_bw/TMath::Sqrt(SBRcounts) << endl;
          aveSBRcount+=signal_counts/TMath::Sqrt(SBRcounts);
          aveSBRbw+=signal_counts_bw/TMath::Sqrt(SBRcounts);
        }

	TF1 *f_cos_count = new TF1("f_cos_count",SpinDensity,0.0,1.0,2);
	f_cos_count->SetParameter(0,0.33);
	f_cos_count->SetParameter(1,1000);
	h_mCounts[KEY_Count]->Fit(f_cos_count,"NQMRI");
	ParSpin_Count[KEY_Count].clear();
	ParSpin_Count[KEY_Count].push_back(static_cast<float>(f_cos_count->GetParameter(0)));
	ParSpin_Count[KEY_Count].push_back(static_cast<float>(f_cos_count->GetParameter(1)));
	h_mRho00[KEY_Rho00_Count]->SetBinContent(h_mRho00[KEY_Rho00_Count]->FindBin(pt_mean),f_cos_count->GetParameter(0));
	h_mRho00[KEY_Rho00_Count]->SetBinError(h_mRho00[KEY_Rho00_Count]->FindBin(pt_mean),f_cos_count->GetParError(0));

	TF1 *f_cos_bw = new TF1("f_cos_bw",SpinDensity,0.0,1.0,2);
	f_cos_bw->SetParameter(0,0.33);
	f_cos_bw->SetParameter(1,1000);
	h_mCounts[KEY_BW]->Fit(f_cos_bw,"NQMI");
	ParSpin_BW[KEY_BW].clear();
	ParSpin_BW[KEY_BW].push_back(static_cast<float>(f_cos_bw->GetParameter(0)));
	ParSpin_BW[KEY_BW].push_back(static_cast<float>(f_cos_bw->GetParameter(1)));
	h_mRho00[KEY_Rho00_BW]->SetBinContent(h_mRho00[KEY_Rho00_BW]->FindBin(pt_mean),f_cos_bw->GetParameter(0));
	h_mRho00[KEY_Rho00_BW]->SetBinError(h_mRho00[KEY_Rho00_BW]->FindBin(pt_mean),f_cos_bw->GetParError(0));
      }
      if(i_cent == 9)
      {
        cout << "Average COUNTS Significance = " << aveSBRcount/totalSBR << endl;
        //cout << "Low COUNTS Significance = " << lowSigCount << endl;
        cout << "Average INTEGRAL Significance = " << aveSBRbw/totalSBR << endl;
        //cout << "Low INTEGRAL Significance = " << lowSigBW << endl;
      }
  }

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",ptQA,9,vmsa::mPID[pid].c_str());
  TCanvas *c_mMass_psi = new TCanvas("c_mMass_psi","c_mMass_psi",10,10,800,800);
  string Title_CosThetaStar[3] = {"2/7 < cos(#theta*) < 3/7","3/7 < cos(#theta*) < 4/7","4/7 < cos(#theta*) < 5/7"};
  c_mMass_psi->Divide(2,2);
  for(int i_theta = 2; i_theta < 5; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi->cd(i_theta-1)->SetLeftMargin(0.15);
    c_mMass_psi->cd(i_theta-1)->SetBottomMargin(0.15);
    c_mMass_psi->cd(i_theta-1)->SetTicks(1,1);
    c_mMass_psi->cd(i_theta-1)->SetGrid(0,0);
    c_mMass_psi->cd(i_theta-1);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",ptQA,9,i_theta,vmsa::mPID[pid].c_str());
    // h_mMass[KEY_QA]->SetTitle(Title_CosThetaStar[i_theta-2].c_str());
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
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]/7.0);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    h_mMass[KEY_QA]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    float x1 = ParBW[KEY_theta_QA][0] - nSigVec*ParBW[KEY_theta_QA][1];
    float x2 = ParBW[KEY_theta_QA][0] + nSigVec*ParBW[KEY_theta_QA][1];
    float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(ParBW[KEY_theta_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    PlotLine(0.6,1.2,0.0,0.0,1,2,2);
  }

  c_mMass_psi->cd(4);
  c_mMass_psi->cd(4)->SetLeftMargin(0.15);
  c_mMass_psi->cd(4)->SetBottomMargin(0.15);
  c_mMass_psi->cd(4)->SetTicks(1,1);
  c_mMass_psi->cd(4)->SetGrid(0,0);
  string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",ptQA,9,vmsa::mPID[pid].c_str());
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
  TF1 *f_cos_count_QA = new TF1("f_cos_count_QA",SpinDensity,0.0,1.0,2);
  f_cos_count_QA->FixParameter(0,ParSpin_Count[KEY_Count][0]);
  f_cos_count_QA->FixParameter(1,ParSpin_Count[KEY_Count][1]);
  f_cos_count_QA->SetLineStyle(2);
  f_cos_count_QA->SetLineColor(4);
  f_cos_count_QA->DrawCopy("l same");

  string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",ptQA,9,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(1.2);
  h_mCounts[KEY_BW]->DrawCopy("pE same");
  TF1 *f_cos_bw_QA = new TF1("f_cos_bw_QA",SpinDensity,0.0,1.0,2);
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
  c_mMass_psi->SaveAs("./figures/KStar/phi_SpinAlighment.png");
#endif
#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  //string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  TCanvas *c_mMass_psi_4 = new TCanvas("c_mMass_psi_4","c_mMass_psi_4",10,10,900,900);
  c_mMass_psi_4->Divide(3,3);
  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi_4->cd(i_theta+1)->SetLeftMargin(0.15);
    c_mMass_psi_4->cd(i_theta+1)->SetBottomMargin(0.15);
    c_mMass_psi_4->cd(i_theta+1)->SetTicks(1,1);
    c_mMass_psi_4->cd(i_theta+1)->SetGrid(0,0);
    c_mMass_psi_4->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",ptQA,9,i_theta,vmsa::mPID[pid].c_str());
    h_mMassBG[KEY_QA]->SetTitle(Form("%d/7 < cos(#theta*) < %d/7",i_theta,i_theta+1));
    h_mMassBG[KEY_QA]->SetMarkerColor(1);
    h_mMassBG[KEY_QA]->SetMarkerStyle(24);
    h_mMassBG[KEY_QA]->SetMarkerSize(0.8);
    h_mMassBG[KEY_QA]->DrawCopy("pE");

    TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
    TF1 *f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
    f_bw->FixParameter(0,ParFit_theta[KEY_theta_QA][0]);
    f_bw->FixParameter(1,ParFit_theta[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParFit_theta[KEY_theta_QA][2]/7.0);
    f_bw->SetParameter(3,ParFit_theta[KEY_theta_QA][3]/7.0);
    f_bw->SetParameter(4,ParFit_theta[KEY_theta_QA][4]/7.0);
    f_bw->SetParameter(5,ParFit_theta[KEY_theta_QA][5]/7.0);
    f_bw->SetParameter(6,ParFit_theta[KEY_theta_QA][6]/7.0);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    h_mMassBG[KEY_QA]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");
    f_bg->FixParameter(0,f_bw->GetParameter(3));
    f_bg->FixParameter(1,f_bw->GetParameter(4));
    f_bg->FixParameter(2,f_bw->GetParameter(5));
    f_bg->FixParameter(3,f_bw->GetParameter(6));
    f_bg->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    f_bg->SetLineColor(kBlue);
    f_bg->DrawCopy("l same");

    float x1 = ParBW[KEY_theta_QA][0] - nSigVec*ParBW[KEY_theta_QA][1];
    float x2 = ParBW[KEY_theta_QA][0] + nSigVec*ParBW[KEY_theta_QA][1];
    float y = h_mMassBG[KEY_QA]->GetBinContent(h_mMassBG[KEY_QA]->FindBin(ParBW[KEY_theta_QA][0]));
    float yl = h_mMassBG[KEY_QA]->GetMinimum();
    PlotLine(x1,x1,yl,y,4,2,2);
    PlotLine(x2,x2,yl,y,4,2,2);
    PlotLine(0.6,1.2,0.0,0.0,1,2,2);
  }

  c_mMass_psi_4->cd(8);
  h_mMass_BGtotal[KEY_theta_QA]->SetTitle("Integrated Yields");
  h_mMass_BGtotal[KEY_theta_QA]->SetMarkerColor(1);
  h_mMass_BGtotal[KEY_theta_QA]->SetMarkerStyle(24);
  h_mMass_BGtotal[KEY_theta_QA]->SetMarkerSize(0.8);
  h_mMass_BGtotal[KEY_theta_QA]->DrawCopy("pE");
  TF1 *f_bw_QA2 = new TF1("f_bw_QA2",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
  TF1 *f_bg_QA2 = new TF1("f_bg_QA2",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
  f_bw_QA2->FixParameter(0,ParFit_theta[KEY_theta_QA][0]);
  f_bw_QA2->FixParameter(1,ParFit_theta[KEY_theta_QA][1]);
  f_bw_QA2->FixParameter(2,ParFit_theta[KEY_theta_QA][2]);
  f_bw_QA2->FixParameter(3,ParFit_theta[KEY_theta_QA][3]);
  f_bw_QA2->FixParameter(4,ParFit_theta[KEY_theta_QA][4]);
  f_bw_QA2->FixParameter(5,ParFit_theta[KEY_theta_QA][5]);
  f_bw_QA2->FixParameter(6,ParFit_theta[KEY_theta_QA][6]);
  f_bw_QA2->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
  f_bw_QA2->SetLineColor(2);
  f_bw_QA2->Draw("l same");
  string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][ptQA],vmsa::pt_up[energy][ptQA]);
  plotTopLegend((char*)pT_range.c_str(),0.15,0.7,0.08,1,0.0,42,1);
  PlotLine(0.6,1.2,0.0,0.0,1,2,2);
  f_bg_QA2->FixParameter(0,f_bw_QA2->GetParameter(3));
  f_bg_QA2->FixParameter(1,f_bw_QA2->GetParameter(4));
  f_bg_QA2->FixParameter(2,f_bw_QA2->GetParameter(5));
  f_bg_QA2->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
  f_bg_QA2->SetLineColor(kBlue);
  f_bg_QA2->DrawCopy("l same");


  c_mMass_psi_4->cd(9);
  //string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_Count]->SetStats(0);
  h_mCounts[KEY_Count]->SetTitle("20-60%");
  h_mCounts[KEY_Count]->SetTitleSize(0.08);
  h_mCounts[KEY_Count]->SetLineColor(4);
  h_mCounts[KEY_Count]->SetMarkerColor(4);
  h_mCounts[KEY_Count]->SetMarkerStyle(24);
  h_mCounts[KEY_Count]->SetMarkerSize(0.8);
  h_mCounts[KEY_Count]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCounts[KEY_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[KEY_Count]->GetXaxis()->CenterTitle();
  h_mCounts[KEY_Count]->DrawCopy("pE");
  TF1 *f_cos_count_QA_3 = new TF1("f_cos_count_QA_3",SpinDensity,0.0,1.0,2);
  f_cos_count_QA_3->FixParameter(0,ParSpin_Count[KEY_Count][0]);
  f_cos_count_QA_3->FixParameter(1,ParSpin_Count[KEY_Count][1]);
  f_cos_count_QA_3->SetLineStyle(2);
  f_cos_count_QA_3->SetLineColor(4);
  f_cos_count_QA_3->DrawCopy("l same");

  //string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(0.8);
  h_mCounts[KEY_BW]->DrawCopy("pE same");
  TF1 *f_cos_bw_QA_3 = new TF1("f_cos_bw_QA_3",SpinDensity,0.0,1.0,2);
  f_cos_bw_QA_3->FixParameter(0,ParSpin_BW[KEY_BW][0]);
  f_cos_bw_QA_3->FixParameter(1,ParSpin_BW[KEY_BW][1]);
  f_cos_bw_QA_3->SetLineStyle(2);
  f_cos_bw_QA_3->SetLineColor(2);
  f_cos_bw_QA_3->DrawCopy("l same");

  TLegend *leg_temp_3 = new TLegend(0.2,0.6,0.5,0.8);
  leg_temp_3->SetFillColor(10);
  leg_temp_3->SetBorderSize(0.0);
  leg_temp_3->AddEntry(h_mCounts[KEY_Count],"bin counting","p");
  leg_temp_3->AddEntry(h_mCounts[KEY_BW],"breit wigner","p");
  leg_temp_3->Draw("same");
  c_mMass_psi_4->SaveAs("./figures/KStar/phi_SpinAlighment_noBGsub.pdf");
#endif


#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  //string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  /*TCanvas *c_mMass_psi_2 = new TCanvas("c_mMass_psi_2","c_mMass_psi_2",10,10,900,900);
  c_mMass_psi_2->Divide(3,3);
  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi_2->cd(i_theta+1)->SetLeftMargin(0.15);
    c_mMass_psi_2->cd(i_theta+1)->SetBottomMargin(0.15);
    c_mMass_psi_2->cd(i_theta+1)->SetTicks(1,1);
    c_mMass_psi_2->cd(i_theta+1)->SetGrid(0,0);
    c_mMass_psi_2->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",ptQA,9,i_theta,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->SetTitle(Form("%d/7 < cos(#theta*) < %d/7",i_theta,i_theta+1));
    h_mMass[KEY_QA]->SetMarkerColor(1);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(0.8);
    h_mMass[KEY_QA]->DrawCopy("pE");

    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
    f_bw->FixParameter(0,ParBW[KEY_theta_QA][0]);
    f_bw->FixParameter(1,ParBW[KEY_theta_QA][1]);
    f_bw->SetParameter(2,ParBW[KEY_theta_QA][2]/7.0);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    h_mMass[KEY_QA]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    float x1 = ParBW[KEY_theta_QA][0] - nSigVec*ParBW[KEY_theta_QA][1];
    float x2 = ParBW[KEY_theta_QA][0] + nSigVec*ParBW[KEY_theta_QA][1];
    float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(ParBW[KEY_theta_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    PlotLine(0.6,1.2,0.0,0.0,1,2,2);
  }

  c_mMass_psi_2->cd(8);
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
  //string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][2],vmsa::pt_up[energy][2]);
  plotTopLegend((char*)pT_range.c_str(),0.15,0.7,0.08,1,0.0,42,1);
  PlotLine(0.6,1.2,0.0,0.0,1,2,2);


  c_mMass_psi_2->cd(9);
  //string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_Count]->SetStats(0);
  h_mCounts[KEY_Count]->SetTitle("20-60%");
  h_mCounts[KEY_Count]->SetTitleSize(0.08);
  h_mCounts[KEY_Count]->SetLineColor(4);
  h_mCounts[KEY_Count]->SetMarkerColor(4);
  h_mCounts[KEY_Count]->SetMarkerStyle(24);
  h_mCounts[KEY_Count]->SetMarkerSize(0.8);
  h_mCounts[KEY_Count]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCounts[KEY_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[KEY_Count]->GetXaxis()->CenterTitle();
  h_mCounts[KEY_Count]->DrawCopy("pE");
  TF1 *f_cos_count_QA_2 = new TF1("f_cos_count_QA_2",SpinDensity,0.0,1.0,2);
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
  TF1 *f_cos_bw_QA_2 = new TF1("f_cos_bw_QA_2",SpinDensity,0.0,1.0,2);
  f_cos_bw_QA_2->FixParameter(0,ParSpin_BW[KEY_BW][0]);
  f_cos_bw_QA_2->FixParameter(1,ParSpin_BW[KEY_BW][1]);
  f_cos_bw_QA_2->SetLineStyle(2);
  f_cos_bw_QA_2->SetLineColor(2);
  f_cos_bw_QA_2->DrawCopy("l same");

  TLegend *leg_temp_2 = new TLegend(0.2,0.6,0.5,0.8);
  leg_temp_2->SetFillColor(10);
  leg_temp_2->SetBorderSize(0.0);
  leg_temp_2->AddEntry(h_mCounts[KEY_Count],"bin counting","p");
  leg_temp_2->AddEntry(h_mCounts[KEY_BW],"breit wigner","p");
  leg_temp_2->Draw("same");
  c_mMass_psi_2->SaveAs("./figures/KStar/phi_SpinAlighment.pdf");
*/
#endif

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  //string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  TH1FMap h_mMassR;
  TH1FMap h_mMassR_total;

  TCanvas *c_mMass_psi_3 = new TCanvas("c_mMass_psi_3","c_mMass_psi_3",10,10,900,900);
  c_mMass_psi_3->Divide(3,3);
  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
  {
    c_mMass_psi_3->cd(i_theta+1)->SetLeftMargin(0.15);
    c_mMass_psi_3->cd(i_theta+1)->SetBottomMargin(0.15);
    c_mMass_psi_3->cd(i_theta+1)->SetTicks(1,1);
    c_mMass_psi_3->cd(i_theta+1)->SetGrid(0,0);
    c_mMass_psi_3->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_%s_SM",ptQA,9,i_theta,vmsa::mPID[pid].c_str());
    h_mMassR[KEY_QA] = (TH1F*)h_mMass[KEY_QA]->Clone();
    h_mMassR[KEY_QA]->SetTitle(Form("%d/7 < cos(#theta*) < %d/7",i_theta,i_theta+1));
    h_mMassR[KEY_QA]->SetMarkerColor(1);
    h_mMassR[KEY_QA]->SetMarkerStyle(24);
    h_mMassR[KEY_QA]->SetMarkerSize(0.8);

    TF1 *f_bwR = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
    f_bwR->FixParameter(0,ParBW[KEY_theta_QA][0]);
    f_bwR->FixParameter(1,ParBW[KEY_theta_QA][1]);
    f_bwR->SetParameter(2,ParBW[KEY_theta_QA][2]/7.0);
    f_bwR->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    h_mMass[KEY_QA]->Fit(f_bwR,"NMQR");
    f_bwR->SetLineColor(2);
    //f_bwR->DrawCopy("l same");
    h_mMassR[KEY_QA]->Divide(f_bwR);
    h_mMassR[KEY_QA]->GetYaxis()->SetRangeUser(0.9,1.1);
    h_mMassR[KEY_QA]->GetXaxis()->SetRangeUser(1.01,1.03);
    h_mMassR[KEY_QA]->DrawCopy("pE");

    float x1 = ParBW[KEY_theta_QA][0] - nSigVec*ParBW[KEY_theta_QA][1];
    float x2 = ParBW[KEY_theta_QA][0] + nSigVec*ParBW[KEY_theta_QA][1];
    //float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(ParBW[KEY_theta_QA][0]));
    //PlotLine(x1,x1,0,1.0,4,2,2);
    //PlotLine(x2,x2,0,1.0,4,2,2);
    PlotLine(0.6,1.2,1.0,1.0,1,2,2);
  }

  c_mMass_psi_3->cd(8);
  h_mMassR_total[KEY_theta_QA] = (TH1F*)h_mMass_total[KEY_theta_QA]->Clone();
  h_mMassR_total[KEY_theta_QA]->SetTitle("Integrated Yield Ratios");
  h_mMassR_total[KEY_theta_QA]->SetMarkerColor(1);
  h_mMassR_total[KEY_theta_QA]->SetMarkerStyle(24);
  h_mMassR_total[KEY_theta_QA]->SetMarkerSize(0.8);
  TF1 *f_bw_QAR = new TF1("f_bw_QAR",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);;
  f_bw_QAR->FixParameter(0,ParBW[KEY_theta_QA][0]);
  f_bw_QAR->FixParameter(1,ParBW[KEY_theta_QA][1]);
  f_bw_QAR->FixParameter(2,ParBW[KEY_theta_QA][2]);
  f_bw_QAR->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
  f_bw_QAR->SetLineColor(2);
  //f_bw_QAR->Draw("l same");
  h_mMassR_total[KEY_theta_QA]->Divide(f_bw_QAR);
  h_mMassR_total[KEY_theta_QA]->GetYaxis()->SetRangeUser(0.9,1.1);
  h_mMassR_total[KEY_theta_QA]->GetXaxis()->SetRangeUser(1.01,1.03);
  h_mMassR_total[KEY_theta_QA]->DrawCopy("pE");
  string pT_rangeR = Form("[%.2f,%.2f]",vmsa::pt_low[energy][2],vmsa::pt_up[energy][2]);
  plotTopLegend((char*)pT_rangeR.c_str(),0.15,0.7,0.08,1,0.0,42,1);
  PlotLine(0.6,1.2,1.0,1.0,1,2,2);


  c_mMass_psi_3->cd(9);
  //string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_Count]->SetStats(0);
  h_mCounts[KEY_Count]->SetTitle("20-60%");
  h_mCounts[KEY_Count]->SetTitleSize(0.08);
  h_mCounts[KEY_Count]->SetLineColor(4);
  h_mCounts[KEY_Count]->SetMarkerColor(4);
  h_mCounts[KEY_Count]->SetMarkerStyle(24);
  h_mCounts[KEY_Count]->SetMarkerSize(0.8);
  h_mCounts[KEY_Count]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCounts[KEY_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[KEY_Count]->GetXaxis()->CenterTitle();
  h_mCounts[KEY_Count]->DrawCopy("pE");
  TF1 *f_cos_count_QA_2R = new TF1("f_cos_count_QA_2R",SpinDensity,0.0,1.0,2);
  f_cos_count_QA_2R->FixParameter(0,ParSpin_Count[KEY_Count][0]);
  f_cos_count_QA_2R->FixParameter(1,ParSpin_Count[KEY_Count][1]);
  f_cos_count_QA_2R->SetLineStyle(2);
  f_cos_count_QA_2R->SetLineColor(4);
  f_cos_count_QA_2R->DrawCopy("l same");

  //string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",2,vmsa::Cent_start,vmsa::mPID[pid].c_str());
  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(0.8);
  h_mCounts[KEY_BW]->DrawCopy("pE same");
  TF1 *f_cos_bw_QA_2R = new TF1("f_cos_bw_QA_2R",SpinDensity,0.0,1.0,2);
  f_cos_bw_QA_2R->FixParameter(0,ParSpin_BW[KEY_BW][0]);
  f_cos_bw_QA_2R->FixParameter(1,ParSpin_BW[KEY_BW][1]);
  f_cos_bw_QA_2R->SetLineStyle(2);
  f_cos_bw_QA_2R->SetLineColor(2);
  f_cos_bw_QA_2R->DrawCopy("l same");

  TLegend *leg_temp_2R = new TLegend(0.5,0.6,0.8,0.8);
  leg_temp_2R->SetFillColor(10);
  leg_temp_2R->SetBorderSize(0.0);
  leg_temp_2R->AddEntry(h_mCounts[KEY_Count],"bin counting","p");
  leg_temp_2R->AddEntry(h_mCounts[KEY_BW],"breit wigner","p");
  leg_temp_2R->Draw("same");
  c_mMass_psi_3->SaveAs("./figures/KStar/phi_SpinYields_Ratios.pdf");
#endif

  // set final pt and rho00 to one TGraphAsymmErrors
  double weightedMeanC = 0.0;
  double weightedErrorC = 0.0;
  double weightedMeanI = 0.0;
  double weightedErrorI = 0.0;
  
  TGraMap g_mRho00;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    //{
      string KEY_Rho00_Count = Form("Rho00_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      g_mRho00[KEY_Rho00_Count] = new TGraphAsymmErrors();
      string KEY_Rho00_BW = Form("Rho00_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
      g_mRho00[KEY_Rho00_BW] = new TGraphAsymmErrors();
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++)
      {
	float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

	float count_content = h_mRho00[KEY_Rho00_Count]->GetBinContent(h_mRho00[KEY_Rho00_Count]->FindBin(pt_mean)); // bin counting
	float count_error   = h_mRho00[KEY_Rho00_Count]->GetBinError(h_mRho00[KEY_Rho00_Count]->FindBin(pt_mean));
	g_mRho00[KEY_Rho00_Count]->SetPoint(i_pt,pt_mean+0.05,count_content);
	g_mRho00[KEY_Rho00_Count]->SetPointError(i_pt,0.0,0.0,count_error,count_error);
	string Name_count = Form("Rho00_%s_Centrality_%d_Count",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
	g_mRho00[KEY_Rho00_Count]->SetName(Name_count.c_str());
	g_mRho00[KEY_Rho00_Count]->SetMarkerStyle(24);
	g_mRho00[KEY_Rho00_Count]->SetMarkerColor(4);

	float bw_content = h_mRho00[KEY_Rho00_BW]->GetBinContent(h_mRho00[KEY_Rho00_BW]->FindBin(pt_mean)); // breit wigner fits
	float bw_error   = h_mRho00[KEY_Rho00_BW]->GetBinError(h_mRho00[KEY_Rho00_BW]->FindBin(pt_mean));
	g_mRho00[KEY_Rho00_BW]->SetPoint(i_pt,pt_mean,bw_content);
	g_mRho00[KEY_Rho00_BW]->SetPointError(i_pt,0.0,0.0,bw_error,bw_error);
	string Name_bw = Form("Rho00_%s_Centrality_%d_BW",vmsa::mPID[pid].c_str(),i_cent); // gaussian fits
	g_mRho00[KEY_Rho00_BW]->SetName(Name_bw.c_str());
	g_mRho00[KEY_Rho00_BW]->SetMarkerStyle(24);
	g_mRho00[KEY_Rho00_BW]->SetMarkerColor(2);
        if(i_cent == 9)
        {
          weightedMeanC  += (count_content/(count_error*count_error));
          weightedErrorC += (1.0/(count_error*count_error));
          weightedMeanI  += (bw_content/(bw_error*bw_error));
          weightedErrorI += (1.0/(bw_error*bw_error));
        
        }
      }
    //}
  }

  std::cout << "Weighted Mean (Count)    = " << weightedMeanC/weightedErrorC << " +/- " << TMath::Sqrt(1.0/weightedErrorC) << std::endl;
  std::cout << "Weighted Mean (Integral) = " << weightedMeanI/weightedErrorI << " +/- " << TMath::Sqrt(1.0/weightedErrorI) << std::endl;
  
  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  // QA Plots for rho00 vs. pt of bin counting and breit wigner integrating
  TCanvas *c_Rho00_pT = new TCanvas("c_Rho00_pT","c_Rho00_pT",10,10,800,800);
  c_Rho00_pT->cd();
  c_Rho00_pT->cd()->SetLeftMargin(0.15);
  c_Rho00_pT->cd()->SetBottomMargin(0.15);
  c_Rho00_pT->cd()->SetTicks(1,1);
  c_Rho00_pT->cd()->SetGrid(0,0);
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

  h_play->GetYaxis()->SetRangeUser(0.2,0.6);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("#rho_{00}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  string KEY_Rho00_Count_QA = Form("Rho00_Count_Centrality_%d_%s",9,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_Count_QA] ,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.447,0.03,1,0.0,42,0);

  string KEY_Rho00_BW_QA = Form("Rho00_BW_Centrality_%d_%s",9,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_BW_QA] ,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  //string leg_besi = "BES-I (counting)";
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19,25,1,1.1);
  //Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,25,1,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  PlotLine(3.0,3.5,0.45,0.45,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),3.6,0.447,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),1.4,0.27,0.04,1,0.0,42,0);
  string leg_centrality = "20%-60%";
  plotTopLegend((char*)leg_centrality.c_str(),2.0,0.24,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/KStar/rho00_pT_%s.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_Rho00_pT->SaveAs(figure_name.c_str());

  
  string OutPutFile = "../output/KStar/RawPhiPt/RawPhiPtSys.root"; 
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_play->Write();
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    //for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    //{
      string KEY_Rho00_Count = Form("Rho00_Count_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // bin counting
      g_mRho00[KEY_Rho00_Count]->Write();
      string KEY_Rho00_BW = Form("Rho00_BW_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // breit wigner fits
      g_mRho00[KEY_Rho00_BW]->Write();
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
      {
	string KEY_Count = Form("Count_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_Count]->Write();
	string KEY_BW = Form("BW_pt_%d_Centrality_%d_2nd_%s_SM",i_pt,i_cent,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_BW]->Write();
      }
    //}
  }
  File_OutPut->Close();
}
