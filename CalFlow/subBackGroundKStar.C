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
//#ifdef MAKECINT
//#pragma link C++ class std::map<std::string,TH1F*>+;
//#endif

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

#ifndef _SaveQA_
#define _SaveQA_  0
#endif

using namespace std;

void subBackGroundKStar(int energy = 4, int pid = 2, int year = 0, int dcaQA = 0, int nsigQA = 0, int nhitQA = 0)
{
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //////// RESOLUTION /////////
  string inputfile = Form("../data/file_%s_Resolution_20220330.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Res = TFile::Open(inputfile.c_str());

  float mTpcSubRes2Val[9];
  float mTpcSubRes2Err[9];
  float mTpcFullRes2Val[9];
  float mTpcFullRes2Err[9];

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
  //////// RESOLUTION /////////

  //string InPutFile_SE = "../data/Yields_KStar_SE_19GeV_20220708.root"; // TPC || TOF
  string InPutFile_SE = "../data/Yields_KStar_SE_19GeV_20221116_TPConly.root"; // TPC only
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  //string InPutFile_ME = "../data/Yields_KStar_ME_19GeV_20220708.root"; // TPC || TOF
  string InPutFile_ME = "../data/Yields_KStar_ME_19GeV_20221116_TPConly.root"; // TPC only
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // Centrality Definitions: cent 0 - 8 are the standard 9 centrality bins
  // cent = 9 : 0-80%, cent = 10 : 0-10%, cent = 11 : 10-40%, cent = 12 : 40-80%

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(Int_t i_pt = 0; i_pt < vmsa::pt_totalKS; i_pt++) // pt bin 
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // centrality bin
    {
      for(Int_t i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // phi-psi bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
            if( i_dca != 0 && i_sig != 0 ) continue;
	    for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
	    {
              if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
	      string KEY_InPutSE = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_SE",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str());
              h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 
              //h_mInPut_SE[KEY_InPutSE]->Rebin(2);
	      string KEY_InPutME = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_ME",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str());
	      h_mInPut_ME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 
              //h_mInPut_ME[KEY_InPutME]->Rebin(2);

	      for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	      {
	        string KEY_SE = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SE",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());
	        string KEY_ME = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_ME",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        h_mMass_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_ME.c_str());
	        string KEY_SM = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	        if(i_norm < 2)
	        {
	          int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_norm]);
	          int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_norm]);

	          float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
	          float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);

	          if(Inte_ME != 0.0) h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
	          h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
	        }
	        else
	        {
	          float Inte_SE = 0.0;
	          float Inte_ME = 0.0;

	          for(int i_inte = 0; i_inte < 2; ++i_inte)
	          {
	            int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_inte]);
	            int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_inte]);
	            Inte_SE += h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
	            Inte_ME += h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	          }

	          h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
	          h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
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

  for(int icent = 0; icent < 13; icent++)
  {

    TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
    c_peak->cd();
    c_peak->cd()->SetLeftMargin(0.15);
    c_peak->cd()->SetBottomMargin(0.15);
    c_peak->cd()->SetTicks(1,1);
    c_peak->cd()->SetGrid(0,0);
    
    string outputname = Form("./figures/KStar/19GeV/allPtBins_%d.pdf",icent);
    string output_start = Form("%s[",outputname.c_str());
    c_peak->Print(output_start.c_str());

     
    for(Int_t i_pt = 0; i_pt < vmsa::pt_totalKS; i_pt++) // pt bin 
    {
      string KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SE",i_pt,icent,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
      h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

      string KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_ME",i_pt,icent,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
      h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
      h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
      h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");


      auto legend = new TLegend(0.5,0.2,0.6,0.4);
      //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
      legend->AddEntry(h_mMass_ME[KEY_ME_QA],"ME bg","f"); 
      legend->AddEntry(h_mMass_SE[KEY_SE_QA],"SE sig+bg","l"); 
      legend->Draw("same");

      if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
      {
        PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
        PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
      }
      if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
      {
        PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
        PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
      }

      c_peak->Update();
      c_peak->Print(outputname.c_str());


      string KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_pt,icent,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_SM[KEY_SM_QA]->GetYaxis()->SetRangeUser(-0.2*h_mMass_SM[KEY_SM_QA]->GetMaximum(),1.1*h_mMass_SM[KEY_SM_QA]->GetMaximum());
      h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
      h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
      h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
      h_mMass_SM[KEY_SM_QA]->DrawCopy("h");

      c_peak->Update();
      c_peak->Print(outputname.c_str());
    }
    string output_stop = Form("%s]",outputname.c_str());
    c_peak->Print(output_stop.c_str()); // close pdf file
  }

  
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
    string KEY_SE_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SE",i_pt,vmsa::Cent_start,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
    h_mMass_SE[KEY_SE_QA]->DrawCopy();

    string KEY_ME_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_ME",i_pt,vmsa::Cent_start,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_pt,vmsa::Cent_start,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStartKS[i_pt],vmsa::ptRawStopKS[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);

    if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
    {
      PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
      PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    }
    if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
    {
      PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
      PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    }
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // Centrality loop
  {
    for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++) // cos(theta*) loop
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
	      for(int pt_bin = vmsa::pt_rebin_firstKSv2[energy]; pt_bin < vmsa::pt_rebin_lastKSv2[energy]; pt_bin++) // pt loop
	      {
	        string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",pt_bin,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        for(int i_pt = vmsa::pt_rebin_startKSv2[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stopKSv2[energy][pt_bin]; i_pt++)
	        {
	          string KEY_SM = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	          //cout << "KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	          if(i_pt == vmsa::pt_rebin_startKSv2[energy][pt_bin])
	          { 
	            cout << "CREATED ==> KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	            h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone(KEY.c_str());
	          }
	          else
	          {
	            h_mMass[KEY]->Add(h_mMass_SM[KEY_SM],1.0);
	            cout << "ADDED TO ==> KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	          }
	        }
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
  for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
  {
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    cout << KEY_QA << endl;
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(1.2);
    h_mMass[KEY_QA]->SetLineColor(kGray+2);
    h_mMass[KEY_QA]->DrawCopy("pE");
    cout << KEY_QA << endl;
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    cout << KEY_QA << endl;
  }
#endif

/*  if(pid == 0 || pid == 2) // Polynomial fit subtraction is only needed for phi meson
  {
    // Poly + Breit Wignar fit to phi integrated InvMass
    TH1FMap h_mMass_theta;
    vecFMap ParFit_theta;
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
	        string KEY_theta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++) // cos(theta*) loop
	        {
	          string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	          if(i_theta == vmsa::PhiPsi_start) h_mMass_theta[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone(KEY_theta.c_str());
	          else h_mMass_theta[KEY_theta]->Add(h_mMass[KEY],1.0);
	        }

	        TF1 *f_sig = new TF1("f_sig",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
	        f_sig->SetParameter(0,vmsa::InvMass[pid]);
	        f_sig->SetParLimits(0,vmsa::InvMass[pid]-1.5*vmsa::Width[pid],vmsa::InvMass[pid]+1.5*vmsa::Width[pid]);
	        f_sig->SetParameter(1,vmsa::Width[pid]);
	        f_sig->SetParLimits(1,0.03,0.065);
	        f_sig->SetParameter(2,1.0);
	        f_sig->SetParameter(3,-1.0);
	        f_sig->SetParameter(4,1.0);
	        f_sig->SetParameter(5,0.1);
	        f_sig->SetParameter(2,h_mMass_theta[KEY_theta]->GetMaximum()/f_sig->GetMaximum());
	        f_sig->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	        //cout << "i_pt = " << i_pt << ", i_cent = " << i_cent << ", i_dca = " << i_dca << ", i_sig = " << i_sig << ", i_norm = " << i_norm << endl;
	        h_mMass_theta[KEY_theta]->Fit(f_sig,"MNR");

	        ParFit_theta[KEY_theta].clear();
	        for(int i_par = 0; i_par < 6; ++i_par)
	        {
	          //cout << "i_par = " << i_par << ", value = " << f_sig->GetParameter(i_par) << endl;;
	          ParFit_theta[KEY_theta].push_back(static_cast<float>(f_sig->GetParameter(i_par)));
	        }
	      }
            }
	  }
	}
      }
    }

#if _PlotQA_
    // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
    TCanvas *c_mMass_theta = new TCanvas("c_mMass_theta","c_mMass_theta",10,10,1400,1400);
    c_mMass_theta->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_theta->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_theta[KEY_theta_QA]->SetMarkerColor(1);
      h_mMass_theta[KEY_theta_QA]->SetMarkerStyle(24);
      h_mMass_theta[KEY_theta_QA]->SetMarkerSize(0.8);
      h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");

      TF1 *f_sig = new TF1("f_sig",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
      for(int i_par = 0; i_par < 6; ++i_par)
      {
	f_sig->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par]);
      }
      f_sig->SetLineColor(2);
      f_sig->SetLineStyle(1);
      f_sig->SetLineWidth(2);
      f_sig->DrawCopy("l same");

      TF1 *f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
      for(int i_par = 0; i_par < 3;++i_par)
      {
	f_bg->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par+3]);
      }
      f_bg->SetLineColor(4);
      f_bg->SetLineStyle(2);
      f_bg->SetLineWidth(4);
      f_bg->DrawCopy("l same");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
#endif

    // Poly+bw fits for phi differential InvMass
    vecFMap ParFit;
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
	        for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++) // cos(theta*) loop
	        {
	          string KEY_theta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	          TF1 *f_sig = new TF1("f_sig",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
	          f_sig->FixParameter(0,ParFit_theta[KEY_theta][0]);
	          f_sig->FixParameter(1,ParFit_theta[KEY_theta][1]);
	          f_sig->SetParameter(2,ParFit_theta[KEY_theta][2]/10.0);
	          f_sig->SetParameter(3,ParFit_theta[KEY_theta][3]/10.0);
	          f_sig->SetParameter(4,ParFit_theta[KEY_theta][4]/10.0);
	          f_sig->SetParameter(5,ParFit_theta[KEY_theta][5]/10.0);
	          f_sig->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);

	          //cout << "i_pt = " << i_pt << ", i_cent = " << i_cent << ", i_theta = " << i_theta << ", i_dca = " << i_dca << ", i_sig = " << i_sig << ", i_norm = " << i_norm << endl;
	          string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	          h_mMass[KEY]->Fit(f_sig,"MNRI");
	          ParFit[KEY].clear();
	          for(int i_par = 0; i_par < 6; ++i_par)
	          {
	            ParFit[KEY].push_back(static_cast<float>(f_sig->GetParameter(i_par)));
	          }

	          TF1 *f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
	          for(int i_par = 0; i_par < 3;++i_par)
	          {
	            f_bg->SetParameter(i_par,ParFit[KEY][i_par+3]);
	          }

	          h_mMass[KEY]->Add(f_bg,-1.0);
	        }
	      }
	  }
	}
      }
    }

#if _PlotQA_
    // QA plots for phi differential InvMass after linear background subtraction
    TCanvas *c_mMass_phi_diff = new TCanvas("c_mMass_phi_diff","c_mMass_phi_diff",10,10,1400,1400);
    c_mMass_phi_diff->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_firstKSv2[energy]; i_pt < vmsa::pt_rebin_lastKSv2[energy]; i_pt++) // pt loop
    {
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_startKSv2[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_QA = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::PhiPsi_start,vmsa::Dca_start,vmsa::nSigKaon_start,0,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("PE");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_lowKSv2[energy][i_pt],vmsa::pt_upKSv2[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
#endif
  }*/
  ////////////////////////////////////////////////////////////////////
  // Yields
  cout << "BEFORE YIELDS" << endl;
  TH1FMap h_mYieldInSE, h_mYieldInME, h_mYield_SE, h_mYield_ME, h_mYield_SM;
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
  { 
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
	{
          if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;

          string KEY_InPutSE = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_SE",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str());
          h_mYieldInSE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone();
          cout << "Loaded --> " << KEY_InPutSE << endl;
          string KEY_InPutME = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_ME",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str());
          h_mYieldInME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone();
          cout << "Loaded --> " << KEY_InPutME << endl;

          for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
          {
            string KEY_SE = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SE",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);  
            h_mYield_SE[KEY_SE] = (TH1F*)h_mYieldInSE[KEY_InPutSE]->Clone(KEY_SE.c_str()); 
            //cout << "Create --> " << KEY_SE << endl;
          
            string KEY_ME = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_ME",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);          
            h_mYield_ME[KEY_ME] = (TH1F*)h_mYieldInME[KEY_InPutME]->Clone(KEY_ME.c_str()); 
            //cout << "Create --> " << KEY_ME << endl;

            string KEY_SM = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);          
            h_mYield_SM[KEY_SM] = (TH1F*)h_mYieldInSE[KEY_InPutSE]->Clone(KEY_SM.c_str()); 
            //cout << "Create --> " << KEY_SM << endl;

            if(i_norm < 2)
            {
              int Norm_bin_start = h_mYieldInSE[KEY_InPutSE]->FindBin(vmsa::Norm_Start[pid][i_norm]);
              int Norm_bin_stop  = h_mYieldInSE[KEY_InPutSE]->FindBin(vmsa::Norm_Stop[pid][i_norm]);

              float Inte_SE = h_mYieldInSE[KEY_InPutSE]->Integral(Norm_bin_start,Norm_bin_stop);
              float Inte_ME = h_mYieldInME[KEY_InPutME]->Integral(Norm_bin_start,Norm_bin_stop);

              if(Inte_ME != 0.0) h_mYield_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
              h_mYield_SM[KEY_SM]->Add(h_mYield_ME[KEY_ME],-1.0);
            }
            else
            {
              float Inte_SE = 0.0;
              float Inte_ME = 0.0;

              for(int i_inte = 0; i_inte < 2; ++i_inte)
              {
                int Norm_bin_start = h_mYieldInSE[KEY_InPutSE]->FindBin(vmsa::Norm_Start[pid][i_inte]);
                int Norm_bin_stop  = h_mYieldInSE[KEY_InPutSE]->FindBin(vmsa::Norm_Stop[pid][i_inte]);
                Inte_SE += h_mYieldInSE[KEY_InPutSE]->Integral(Norm_bin_start,Norm_bin_stop);
                Inte_ME += h_mYieldInME[KEY_InPutME]->Integral(Norm_bin_start,Norm_bin_stop);
              }

              h_mYield_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
              h_mYield_SM[KEY_SM]->Add(h_mYield_ME[KEY_ME],-1.0);
            }
          }
        }
      }
    }
  }


  /*TH1FMap h_mYield; // for QA plot only
  vecFMap ParYield_SM;

  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
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
            string KEY_Yield = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            h_mYield[KEY_Yield] = (TH1F*)h_mYield_SM[KEY_Yield]->Clone();
            TF1 *f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
            for(Int_t i_par = 0; i_par < 6; i_par++)
            {
              f_bw->ReleaseParameter(i_par);
            }
            f_bw->SetParameter(0,vmsa::InvMass[pid]);
            f_bw->SetParLimits(0,vmsa::InvMass[pid]-1.5*vmsa::Width[pid],vmsa::InvMass[pid]+1.5*vmsa::Width[pid]);
            f_bw->SetParameter(1,vmsa::Width[pid]);
            f_bw->SetParameter(2,10000);
            f_bw->SetParameter(3,-6000);
            f_bw->SetParameter(4,0.5);
            f_bw->SetParameter(5,0.1);
            f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
            ParYield_SM[KEY_Yield].clear();
            h_mYield[KEY_Yield]->Fit(f_bw,"QNR");
            for(Int_t n_par = 0; n_par < 6; n_par++)
            {
              ParYield_SM[KEY_Yield].push_back(static_cast<Float_t>(f_bw->GetParameter(n_par)));
            }

            TF1 *f_poly = new TF1("f_poly",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
            f_poly->FixParameter(0,f_bw->GetParameter(3));
            f_poly->FixParameter(1,f_bw->GetParameter(4));
            f_poly->FixParameter(2,f_bw->GetParameter(5));
            h_mYield[KEY_Yield]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
          }
        }
      }
    }
  }*/
  
  vecFMap ParYield_BW;
  vecFMap yields_Gaus, yields_BW;
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
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
            string KEY_Yield = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            TF1 *f_yields_bw = new TF1("f_yields_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
            f_yields_bw->SetParameter(0,vmsa::InvMass[pid]);
            f_yields_bw->SetParLimits(0,vmsa::InvMass[pid]-0.03,vmsa::InvMass[pid]+0.03);
            f_yields_bw->SetParameter(1,vmsa::Width[pid]);
            f_yields_bw->SetParameter(2,h_mYield_SM[KEY_Yield]->GetMaximum());
            f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
            TFitResultPtr result = h_mYield_SM[KEY_Yield]->Fit(f_yields_bw,"MQNRS");
            ParYield_BW[KEY_Yield].clear();
            ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(0)));
            ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(1)));
            ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(2)));
            ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(3)));
            ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(4)));
            ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(5)));

            TF1 *f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
            f_bg->SetParameter(0,f_yields_bw->GetParameter(3));
            f_bg->SetParameter(1,f_yields_bw->GetParameter(4));
            f_bg->SetParameter(2,f_yields_bw->GetParameter(5));
            f_bg->SetParError(0,f_yields_bw->GetParError(3));
            f_bg->SetParError(1,f_yields_bw->GetParError(4));
            f_bg->SetParError(2,f_yields_bw->GetParError(5));

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

            float bin_width = h_mYield_SM[KEY_Yield]->GetBinWidth(1);
            float Inte_start = ParYield_BW[KEY_Yield][0]-vmsa::nSigVecKS*ParYield_BW[KEY_Yield][1]-0.5*bin_width;
            float Inte_stop  = ParYield_BW[KEY_Yield][0]+vmsa::nSigVecKS*ParYield_BW[KEY_Yield][1]+0.5*bin_width;
            float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
            float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width; 

            // counting for guassian
            float counts_gaus = 0.0;
            float errors_gaus = 0.0;
            int bin_start = h_mYield_SM[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]-vmsa::nSigVecKS*ParYield_BW[KEY_Yield][1]);
            int bin_stop  = h_mYield_SM[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]+vmsa::nSigVecKS*ParYield_BW[KEY_Yield][1]);
            for(Int_t i_bin = bin_start; i_bin <= bin_stop; i_bin++)
            {
              counts_gaus += h_mYield_SM[KEY_Yield]->GetBinContent(i_bin);
              errors_gaus += h_mYield_SM[KEY_Yield]->GetBinError(i_bin)*h_mYield_SM[KEY_Yield]->GetBinError(i_bin);
            }
            yields_Gaus[KEY_Yield].clear();
            yields_Gaus[KEY_Yield].push_back(static_cast<Float_t>(counts_gaus-counts_bg));
            yields_Gaus[KEY_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors_gaus+errors_bg*errors_bg)));

            // integrating for breit wigner
            //float bin_width = h_mYield_SM[KEY_Yield]->GetBinWidth(1);
            //float Inte_start = ParYield_BW[KEY_Yield][0]-vmsa::nSigVecKS*ParYield_BW[KEY_Yield][1]-0.5*bin_widith;
            //float Inte_stop  = ParYield_BW[KEY_Yield][0]+vmsa::nSigVecKS*ParYield_BW[KEY_Yield][1]+0.5*bin_width;
            float counts_bw = f_yields_bw->Integral(Inte_start,Inte_stop)/bin_width;
            float errors_bw = f_yields_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
            yields_BW[KEY_Yield].clear();
            yields_BW[KEY_Yield].push_back(static_cast<Float_t>(counts_bw-counts_bg));
            yields_BW[KEY_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg)));
          }
        }
      }
    }
  }

#if _PlotQA_
  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,900);
  c_diff->Divide(3,3);
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    c_diff->cd(i_cent+1);
    c_diff->cd(i_cent+1)->SetLeftMargin(0.15);
    c_diff->cd(i_cent+1)->SetBottomMargin(0.15);
    c_diff->cd(i_cent+1)->SetTicks(1,1);
    c_diff->cd(i_cent+1)->SetGrid(0,0);

    string KEY_QA = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,dcaQA,nsigQA,nhitQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mYield_SM[KEY_QA]->SetTitle("");
    h_mYield_SM[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mYield_SM[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
    h_mYield_SM[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
    h_mYield_SM[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
    h_mYield_SM[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
    h_mYield_SM[KEY_QA]->GetXaxis()->CenterTitle();

    h_mYield_SM[KEY_QA]->GetYaxis()->SetRangeUser(h_mYield_SM[KEY_QA]->GetMinimum(),1.1*h_mYield_SM[KEY_QA]->GetMaximum());
    h_mYield_SM[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
    h_mYield_SM[KEY_QA]->GetYaxis()->SetTitle("Yields");
    h_mYield_SM[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
    h_mYield_SM[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
    h_mYield_SM[KEY_QA]->GetYaxis()->CenterTitle();

    h_mYield_SM[KEY_QA]->SetMarkerStyle(24);
    h_mYield_SM[KEY_QA]->SetMarkerColor(kGray+2);
    h_mYield_SM[KEY_QA]->SetMarkerSize(1.2);
    h_mYield_SM[KEY_QA]->Draw("pE");
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);

    TF1 *f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
    f_bw->SetParameter(0,ParYield_BW[KEY_QA][0]);
    f_bw->SetParameter(1,ParYield_BW[KEY_QA][1]);
    f_bw->SetParameter(2,ParYield_BW[KEY_QA][2]);
    f_bw->SetParameter(3,ParYield_BW[KEY_QA][3]);
    f_bw->SetParameter(4,ParYield_BW[KEY_QA][4]);
    f_bw->SetParameter(5,ParYield_BW[KEY_QA][5]);
    f_bw->SetLineColor(4);
    f_bw->SetLineStyle(2);
    f_bw->SetLineWidth(2);
    f_bw->Draw("l same");
  }
#endif
 

  Float_t yields_total_gaus[4][vmsa::Dca_stop][vmsa::nSigKaon_stop][vmsa::mNHit_stop] = {0.0};
  Float_t yields_total_bw[4][vmsa::Dca_stop][vmsa::nSigKaon_stop][vmsa::mNHit_stop]   = {0.0};
  // calculate final resolution correction factors and correct flow
  for(Int_t i_cent = vmsa::Cent_start; i_cent < 9; i_cent++) // centrality bin loop
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
            string KEY_Yield = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            //cout << KEY_Yield << endl;
            //cout << "yield gaus = " << yields_Gaus[KEY_Yield][0] << endl;
            //cout << "yield bw   = " << yields_BW[KEY_Yield][0]   << endl;

            if(i_cent >= 0 && i_cent <= 8) // calculate resolution and total yields in 0-80% centrality bin
            {
              yields_total_gaus[0][i_dca][i_sig][i_nhit] += yields_Gaus[KEY_Yield][0];
              yields_total_bw[0][i_dca][i_sig][i_nhit]   += yields_BW[KEY_Yield][0];
            }
            if(i_cent == 8 || i_cent == 7) // calculate resolution and total yields in 0-10% centrality bin
            {
              yields_total_gaus[1][i_dca][i_sig][i_nhit] += yields_Gaus[KEY_Yield][0];
              yields_total_bw[1][i_dca][i_sig][i_nhit]   += yields_BW[KEY_Yield][0];
            }
            if(i_cent <= 6 && i_cent >= 4) // calculate resolution and total yields in 10-40% centrality bin
            {
              yields_total_gaus[2][i_dca][i_sig][i_nhit] += yields_Gaus[KEY_Yield][0];
              yields_total_bw[2][i_dca][i_sig][i_nhit]   += yields_BW[KEY_Yield][0];
            }
            if(i_cent <= 3 && i_cent >= 0) // calculate resolution and total yields in 40-80% centrality bin
            {
              yields_total_gaus[3][i_dca][i_sig][i_nhit] += yields_Gaus[KEY_Yield][0];
              yields_total_bw[3][i_dca][i_sig][i_nhit]   += yields_BW[KEY_Yield][0];
            }
          }
        }
      }
    }
  }

  vecFMap mean_res_gaus;
  vecFMap mean_res_bw;
  for(Int_t i_cent = vmsa::Cent_start; i_cent < 4; i_cent++) // centrality bin loop
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
            string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent+9,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            string KEY_WeightBW   = Form("ResWeight_BW_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent+9,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            mean_res_gaus[KEY_WeightGaus].push_back(0.0);
            mean_res_bw[KEY_WeightBW].push_back(0.0);    
          }
        }
      }
    }
  }

  for(Int_t i_cent = vmsa::Cent_start; i_cent < 9; i_cent++) // centrality bin loop
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
            string KEY_Yield = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);

            if(i_cent >= 0 && i_cent <= 8) // calculate resolution and total yields in 0-80% centrality bin
            {
              string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_9_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              string KEY_WeightBW   = Form("ResWeight_BW_Centrality_9_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              mean_res_gaus[KEY_WeightGaus][0] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[0][i_dca][i_sig][i_nhit]);
              mean_res_bw[KEY_WeightBW][0]     += yields_BW[KEY_Yield][0]  /(mTpcSubRes2Val[i_cent]*yields_total_bw[0][i_dca][i_sig][i_nhit]);
              cout << "cent = " << i_cent << ", i_dca = " << i_dca << ", i_sig = " << i_sig << ", i_nhit = " << i_nhit << endl;
              cout << "Gauss Contribution = " << yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[0][i_dca][i_sig][i_nhit]) << endl;
              cout << "BW    Contribution = " << yields_BW[KEY_Yield][0]  /(mTpcSubRes2Val[i_cent]*yields_total_bw[0][i_dca][i_sig][i_nhit]) << endl;      
            }
            if(i_cent == 8 || i_cent == 7) // calculate resolution and total yields in 0-10% centrality bin
            {
              string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_10_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              string KEY_WeightBW   = Form("ResWeight_BW_Centrality_10_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              mean_res_gaus[KEY_WeightGaus][0] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[1][i_dca][i_sig][i_nhit]);
              mean_res_bw[KEY_WeightBW][0]     += yields_BW[KEY_Yield][0]  /(mTpcSubRes2Val[i_cent]*yields_total_bw[1][i_dca][i_sig][i_nhit]);
            }
            if(i_cent <= 6 && i_cent >= 4) // calculate resolution and total yields in 10-40% centrality bin
            {
              string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_11_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              string KEY_WeightBW   = Form("ResWeight_BW_Centrality_11_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              mean_res_gaus[KEY_WeightGaus][0] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[2][i_dca][i_sig][i_nhit]);
              mean_res_bw[KEY_WeightBW][0]     += yields_BW[KEY_Yield][0]  /(mTpcSubRes2Val[i_cent]*yields_total_bw[2][i_dca][i_sig][i_nhit]);
            }
            if(i_cent <= 3 && i_cent >= 0) // calculate resolution and total yields in 40-80% centrality bin
            {
              string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_12_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              string KEY_WeightBW   = Form("ResWeight_BW_Centrality_12_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
              mean_res_gaus[KEY_WeightGaus][0] += yields_Gaus[KEY_Yield][0]/(mTpcSubRes2Val[i_cent]*yields_total_gaus[3][i_dca][i_sig][i_nhit]);
              mean_res_bw[KEY_WeightBW][0]     += yields_BW[KEY_Yield][0]  /(mTpcSubRes2Val[i_cent]*yields_total_bw[3][i_dca][i_sig][i_nhit]);
            }
          }
        }
      }
    }
  }

  for(Int_t i_cent = vmsa::Cent_start; i_cent < 4; i_cent++) // centrality bin loop
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
            string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent+9,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            string KEY_WeightBW   = Form("ResWeight_BW_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent+9,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            //cout << KEY_WeightGaus << endl;
            //cout << "centrality_bin = " << i_cent+9 << ", mean_res_gaus = " << mean_res_gaus[KEY_WeightGaus][0] << endl;
            //cout << KEY_WeightBW << endl;
            //cout << "centrality_bin = " << i_cent+9 << ", mean_res_bw   = " << mean_res_bw[KEY_WeightBW][0]     << endl;
          }
        }
      }
    }
  }



  /// YIELDS ////////////////////////////////////

  string outputfile = Form("../output/AuAu%s/Flow/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
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
	    for(int i_theta = vmsa::PhiPsi_start; i_theta < vmsa::PhiPsi_stop; i_theta++) // cos(theta*) loop
	    {
	      for(int i_norm = vmsa::Norm_start; i_norm < 1; ++i_norm)
	      {
	        string KEY = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
	        h_mMass[KEY]->Write();
                //cout << "Writing --> " << KEY << endl;
	      }
            }
	  }
	}
      }
    }
  }

  for(Int_t i_cent = 0; i_cent < 13; i_cent++) // Centrality loop
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
            string KEY_WeightGaus = Form("ResWeight_Gauss_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            string KEY_WeightBW   = Form("ResWeight_BW_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_SM",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm);
            if(i_cent > 8)
            {
              File_OutPut->WriteObjectAny(&mean_res_gaus[KEY_WeightGaus],"std::vector<float>",KEY_WeightGaus.c_str()); 
              //cout << "Writing --> " << KEY_WeightGaus << endl;
              File_OutPut->WriteObjectAny(&mean_res_bw[KEY_WeightBW],"std::vector<float>",KEY_WeightBW.c_str()); 
              //cout << "Writing --> " << KEY_WeightBW   << endl;
            }
            else
            {
              float w = float(1.0/mTpcSubRes2Val[i_cent]);
              std::vector<float> weight{w};
              File_OutPut->WriteObjectAny(&weight,"std::vector<float>",KEY_WeightGaus.c_str()); 
              //cout << "Writing --> " << KEY_WeightGaus << endl;
              File_OutPut->WriteObjectAny(&weight,"std::vector<float>",KEY_WeightBW.c_str()); 
              //cout << "Writing --> " << KEY_WeightBW   << endl;
            }
          }
        }
      }
    }
  }

  File_OutPut->Close();
}
