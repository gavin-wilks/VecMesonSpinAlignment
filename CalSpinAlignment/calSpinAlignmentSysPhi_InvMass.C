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


double rhoInvMassPoly2(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];
  double pfit = par[6];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*x + par[5]*x*x;
 
  double y = Poly*(1+(BW*pfit)/(1+BW)); 

  return y;
}


void calSpinAlignmentSysPhi_InvMass(int energy = 3, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1")
{
  
  std::string EP[2] = {"","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg_%s_InvMass_Ratios.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  string inputfileInvMass = Form("../output/AuAu%s/%s/InvMassSubBg_%s_InvMass.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/InvMassSubBg_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TFile *File_InPut_InvMass = TFile::Open(inputfileInvMass.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
  TH1FMap h_mMass_InteTheta_InvMass;
  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  //for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  for(int i_norm = 0; i_norm < 1; ++i_norm)
	  {
	    string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)File_InPut->Get(KEY_InteTheta.c_str());
            h_mMass_InteTheta_InvMass[KEY_InteTheta] = (TH1F*)File_InPut_InvMass->Get(KEY_InteTheta.c_str());           
	    TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	    f_bw->SetParameter(0,vmsa::InvMass[pid]);
	    f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.003,vmsa::InvMass[pid]+0.003);
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
	  }
	}
      }
    }
  }

//#if _PlotQA_
//  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,900);
//  c_diff->Divide(3,3);
//  for(int i_theta = 0; i_theta < 9; ++i_theta)
//  {
//    c_diff->cd(i_theta+1);
//    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
//    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
//    c_diff->cd(i_theta+1)->SetTicks(1,1);
//    c_diff->cd(i_theta+1)->SetGrid(0,0);
//    if(i_theta < vmsa::CTS_stop)
//    {
//      string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",vmsa::pt_QA[energy],9,i_theta,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//      h_mMass[KEY_QA]->SetTitle("");
//      h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
//      h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
//      h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
//      h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
//      h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
//      h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
//
//      h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
//      h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
//      h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
//      h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
//      h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
//      h_mMass[KEY_QA]->GetYaxis()->CenterTitle();
//
//      h_mMass[KEY_QA]->SetMarkerStyle(24);
//      h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
//      h_mMass[KEY_QA]->SetMarkerSize(1.2);
//      h_mMass[KEY_QA]->Draw("pE");
//      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
//    }
//    if(i_theta == vmsa::CTS_stop)
//    {
//      string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",vmsa::pt_QA[energy],9,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
//
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Yields");
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();
//
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
//      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
//      TF1 *f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
//      f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][0]);
//      f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][1]);
//      f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][2]);
//      f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA][3]);
//      f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA][4]);
//      f_bw->SetLineColor(2);
//      f_bw->SetLineStyle(2);
//      f_bw->SetLineWidth(2);
//      f_bw->Draw("l same");
//    }
//  }
//#endif

  /*
  // fit theta-differential bin to extract fit parameters for integration
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
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
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	//for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	for(int i_norm = 0; i_norm < 1; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = 1; i_method < 2; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              //cout << KEY_rho << endl;
	      g_mRho[KEY_rho] = new TGraphAsymmErrors();
	      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	      {
		string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
		TF1 *f_bw = new TF1("f_bw",rhoInvMassPoly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
		f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta][0]); // m0
		f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta][1]); // Gamma
		f_bw->FixParameter(2,Par_InteTheta[KEY_InteTheta][2]); // Norm

                /* Other Parameters
                   3 = pol0;
                   4 = pol1;
                   5 = pol2; 
                   6 = pfit;       */

		f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
		h_mMass_InteTheta_InvMass[KEY_InteTheta]->Fit(f_bw,"MQNR");

                //assuming <cos^2(theta*)>^true_B = 1/3
                double rho_B = 1./3.;
                double pfit = f_bw->GetParameter(6); 
                double pfit_err = f_bw->GetParError(6);
                double rhoRaw = rho_B*5./2.*(1+pfit)-1./2.;
                double rhoRaw_err = rho_B*5./2.*pfit_err;
                              
                cout << "pt = [" << vmsa::pt_low[energy][i_pt] << "," << vmsa::pt_up[energy][i_pt] << "] GeV/c,    rho00 raw = " << rhoRaw << " +/- " << rhoRaw_err << endl;

		Par[KEY_InteTheta].clear();
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(0)));
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(1)));
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(2)));
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(3)));
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(4)));
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(5)));
		Par[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(6)));

		float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

		g_mRho[KEY_rho]->SetPoint(i_pt,pt_mean,rhoRaw);
		g_mRho[KEY_rho]->SetPointError(i_pt,0.0,0.0,rhoRaw_err,rhoRaw_err);
	      }
	    }
	  }
	}
      }
    }
  }
  
#if _PlotQA_
  for(int i_cent = 9; i_cent < 10; ++i_cent) // Centrality loop
  {
    string outputname = Form("./figures/%s/%s/pTstudy/allThetaYields%s_%s_Order%d_cent%d_%s_pTdependence_InvMass.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_cent,etamode.c_str());
    string output_start = Form("%s[",outputname.c_str());
    
    //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
    TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,600,900);
    c_diff->Print(output_start.c_str());
    c_diff->Clear();
    c_diff->Divide(2,3);
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
    {
      c_diff->cd(i_pt+1);
      c_diff->cd(i_pt+1)->SetLeftMargin(0.15);
      c_diff->cd(i_pt+1)->SetBottomMargin(0.15);
      c_diff->cd(i_pt+1)->SetTicks(1,1);
      c_diff->cd(i_pt+1)->SetGrid(0,0);
      string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->SetTitle("");
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetXaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->SetTitle(Form("%.2f<p_{T}<%.2f",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]));
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetXaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetXaxis()->CenterTitle();

      if(i_pt == vmsa::pt_rebin_last[energy]-1) h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->SetRangeUser(0.8*h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetMinimum(),1.075*h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetMaximum());
      else h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->SetRangeUser(0.95*h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetMinimum(),1.05*h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetMaximum());
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->SetTitle("Yields");
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->GetYaxis()->CenterTitle();

      h_mMass_InteTheta_InvMass[KEY_InteTheta]->SetMarkerStyle(24);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->SetMarkerColor(kGray+2);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->SetMarkerSize(1.2);
      h_mMass_InteTheta_InvMass[KEY_InteTheta]->Draw("pE");
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);

      TF1 *f_bw = new TF1("f_bw",rhoInvMassPoly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
      f_bw->SetParameter(0,Par[KEY_InteTheta][0]);
      f_bw->SetParameter(1,Par[KEY_InteTheta][1]);
      f_bw->SetParameter(2,Par[KEY_InteTheta][2]);
      f_bw->SetParameter(3,Par[KEY_InteTheta][3]);
      f_bw->SetParameter(4,Par[KEY_InteTheta][4]);
      f_bw->SetParameter(5,Par[KEY_InteTheta][5]);
      f_bw->SetParameter(6,Par[KEY_InteTheta][6]);
      f_bw->SetLineColor(kOrange+7);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->Draw("l same");

      TF1 *f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
      f_bg->SetParameter(0,Par[KEY_InteTheta][3]);
      f_bg->SetParameter(1,Par[KEY_InteTheta][4]);
      f_bg->SetParameter(2,Par[KEY_InteTheta][5]);
      f_bg->SetLineColor(kBlue);
      f_bg->SetLineStyle(2);
      f_bg->SetLineWidth(2);
      f_bg->Draw("l same");

      TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
      leg1->AddEntry(h_mMass_InteTheta_InvMass[KEY_InteTheta],"data","p");
      leg1->AddEntry(f_bw,"sig+res","l");
      leg1->AddEntry(f_bg,"res","l");
      leg1->Draw("same");
    }

    c_diff->Update();
    c_diff->Print(outputname.c_str());
    string output_stop = Form("%s]",outputname.c_str());
    c_diff->Print(output_stop.c_str()); // close pdf file
  }
   // if(!random3D) c_diff->SaveAs("../figures/c_diff_2.pdf");
   // if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2.pdf");
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

 // for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
 // {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < 1; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < 1; ++i_sigma)
	  {
	    for(int i_method = 1; i_method < 2; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",9,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,0,1.1);
	    }
	  }
	}
      }
    }
  //}
  if(!random3D) c_rho->SaveAs("../figures/c_rho.pdf");
  if(random3D)  c_rho->SaveAs("../figures/3DRandom/c_rho.pdf");

  //string outputfile = Form("../output/AuAu%s/%s/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string outputfile = Form("../output/AuAu%s/%s/RawPhiPtSys_%s_InvMass.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/RawPhiPtSys_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
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
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = 1; i_method < 2; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
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
