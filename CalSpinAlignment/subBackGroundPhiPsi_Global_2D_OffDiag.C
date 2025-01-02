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

void subBackGroundPhiPsi_Global_2D_OffDiag(int energy = 4, int pid = 0, int year = 0, string date = "20241021", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int i_phipsi = 0, int frameopt = 0)
{
  std::string frame = "Global";
  if(frameopt == 1) frame = "Helicity";

  std::string EP[2] = {"1st","2nd"};
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //string InPutFile_SE = Form("../data/Yields_Phi_SE_19GeV_20220527.root"); //original eta < 1.0
  string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_%s_2D_Study3.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str(),frame.c_str());
  //if(energy == 3) InPutFile_SE = Form("../data/Yields_Phi_SE_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_%s_2D_Study3_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str(),frame.c_str());
  if(random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_19GeV_20220408.root"); //original eta < 1.0
  string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_%s_2D_Study3.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str(),frame.c_str());
  //if(energy == 3) InPutFile_ME = Form("../data/Yields_Phi_ME_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_%s_2D_Study3_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str(),frame.c_str());
  if(random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;

  for(Int_t i_phi = 0; i_phi < 12; i_phi++) // phi' bin
  {
    for(Int_t i_thetah = 0; i_thetah < 9; i_thetah++) // cos(theta*) bin
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;
          for(Int_t ptbin = 0; ptbin < vmsa::pt_rebin_y_2D; ptbin++)
          {
            for(Int_t i_pt = vmsa::pt_rebin_start_y_2D[energy][ptbin]; i_pt <= vmsa::pt_rebin_stop_y_2D[energy][ptbin]; i_pt++)
            {
              for(Int_t centbin = 0; centbin < vmsa::cent_rebin_total_2D; centbin++)
              {
                for(Int_t i_cent = vmsa::cent_rebin_2D[centbin]; i_cent < vmsa::cent_rebin_2D[centbin+1]; i_cent++)
                {
                  for(Int_t i_phipsiRE = i_phipsi*2; i_phipsiRE <= i_phipsi*2+1; i_phipsiRE++)
                  {
                    string KEY_InPutSE = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_SE",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
                    string KEY_InPutSE_Rebin = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_SE",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
                    cout << KEY_InPutSE << endl;
                    cout << KEY_InPutSE_Rebin << endl;
                    if(i_pt == vmsa::pt_rebin_start_y_2D[energy][ptbin] && i_cent == vmsa::cent_rebin_2D[centbin] && i_phipsiRE == i_phipsi*2) h_mInPut_SE[KEY_InPutSE_Rebin] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 
                    else                                                                                           h_mInPut_SE[KEY_InPutSE_Rebin]->Add((TH1F*)File_SE->Get(KEY_InPutSE.c_str()),1.0); 


                    string KEY_InPutME = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_ME",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
                    string KEY_InPutME_Rebin = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_ME",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
                    cout << KEY_InPutME << endl;
                    cout << KEY_InPutME_Rebin << endl;
                    if(i_pt == vmsa::pt_rebin_start_y_2D[energy][ptbin] && i_cent == vmsa::cent_rebin_2D[centbin] && i_phipsiRE == i_phipsi*2) h_mInPut_ME[KEY_InPutME_Rebin] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 
                    else                                                                                           h_mInPut_ME[KEY_InPutME_Rebin]->Add((TH1F*)File_ME->Get(KEY_InPutME.c_str()),1.0); 
                  }
                } 
              }
            }
          }
          for(Int_t ptbin = 0; ptbin < vmsa::pt_rebin_y_2D; ptbin++)
          {
            for(Int_t centbin = 0; centbin < vmsa::cent_rebin_total_2D; centbin++)
            {
              for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
              {
                string KEY_InPutSE = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_SE",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
                string KEY_SE = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());
                string KEY_InPutME = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_ME",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
                string KEY_ME = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_ME.c_str());
                string KEY_SM = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",ptbin,i_phipsi,centbin,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
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

                  if(Inte_ME != 0.0) h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
                  h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  cout << "Loaded all plots" << endl;

#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
  c_peak->cd();
  c_peak->cd()->SetLeftMargin(0.15);
  c_peak->cd()->SetBottomMargin(0.15);
  c_peak->cd()->SetTicks(1,1);
  c_peak->cd()->SetGrid(0,0);
  string KEY_SE_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",0,i_phipsi,0,vmsa::CTS_start,0,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
  h_mMass_SE[KEY_SE_QA]->SetStats(0);
  h_mMass_SE[KEY_SE_QA]->SetTitle("");
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",0,i_phipsi,0,vmsa::CTS_start,0,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",0,i_phipsi,0,vmsa::CTS_start,0,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

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
  TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
  leg1->AddEntry(h_mMass_SE[KEY_SE_QA],"Same Event","l");
  leg1->AddEntry(h_mMass_ME[KEY_ME_QA],"Mixed Event","f");
  leg1->AddEntry(h_mMass_SM[KEY_SM_QA],"Same-Mixed Event","f");
  leg1->Draw("same");

  c_peak->SaveAs(Form("figures/%s/%s/phipsistudy/%s2D_CloseUpMixedEvent_%s_Order%d_Cent0_y%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_phipsi));
  

  // QA Plots for pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,2800,800);
  c_pT->Divide(7,2);

  string outputname = Form("figures/%s/%s/phipsistudy/%s2D_InvMassDistributions_%s_Order%d_Cent0_y%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_phipsi);
  string output_start = Form("%s[",outputname.c_str());
  
  c_pT->Print(output_start.c_str());
  for(int iphi = 0; iphi < 12; iphi++)
  {
    for(int ithetah = 0; ithetah < 9; ithetah++)
    {
      c_pT->cd(ithetah+1);
      c_pT->cd(ithetah+1)->SetLeftMargin(0.15);
      c_pT->cd(ithetah+1)->SetBottomMargin(0.15);
      c_pT->cd(ithetah+1)->SetTicks(1,1);
      c_pT->cd(ithetah+1)->SetGrid(0,0);
      string KEY_SE_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",0,i_phipsi,0,ithetah,iphi,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
      h_mMass_SE[KEY_SE_QA]->SetTitle(Form("%d#pi/6<#beta<%d#pi/6, %1.1f/4.5<cos(#theta)<%1.1f/4.5",iphi,iphi+1,float(ithetah)-4.5,float(ithetah)-3.5));
      h_mMass_SE[KEY_SE_QA]->DrawCopy();

      string KEY_ME_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",0,i_phipsi,0,ithetah,iphi,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
      h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
      h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
      h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

      string KEY_SM_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",0,i_phipsi,0,ithetah,iphi,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
      h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
      h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
      h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

      string pT_range = Form("[%d#pi/6,%d#pi/6)",i_phipsi,i_phipsi+1);
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
    c_pT->Update();
    c_pT->Print(outputname.c_str());
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_pT->Print(output_stop.c_str()); // close pdf file
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(Int_t i_pt = 0; i_pt < vmsa::pt_rebin_y_2D; i_pt++)
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++) // Centrality loop
    {
      for(int i_phi = vmsa::CTS_start; i_phi < 12/*2*vmsa::CTS_stop*/; i_phi++) // cos(theta*) loop
      {
        for(int i_thetah = vmsa::CTS_start; i_thetah < 9/*2*vmsa::CTS_stop*/; i_thetah++) // cos(theta*) loop
        {
          for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
          {
            for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
            {
              if( i_dca != 0 && i_sig != 0 ) continue;
              for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
              {
                string KEY = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                string KEY_SM = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone(KEY.c_str());
              }
            }
          }
        }
      }
    }
  }
 
#if _PlotQA_
  // QA Plots for pT rebins
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,2800,800);
  c_pT_rebin->Divide(7,2);

  outputname = Form("figures/%s/%s/phipsistudy/%s2D_InvMassBackgroundSubtracted_Rebinned_%s_Order%d_Cent0_y%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_phipsi);
  output_start = Form("%s[",outputname.c_str());
  
  c_pT_rebin->Print(output_start.c_str());

  for(int iphi = 0; iphi < 12/*2*vmsa::CTS_total*/; iphi++)
  {
    for(int ithetah = 0; ithetah < 9/*2*vmsa::CTS_total*/; ithetah++)
    {
      c_pT_rebin->cd(i_phipsi+1);
      c_pT_rebin->cd(i_phipsi+1)->SetLeftMargin(0.15);
      c_pT_rebin->cd(i_phipsi+1)->SetBottomMargin(0.15);
      c_pT_rebin->cd(i_phipsi+1)->SetTicks(1,1);
      c_pT_rebin->cd(i_phipsi+1)->SetGrid(0,0);
      string KEY_QA = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",0,i_phipsi,0,ithetah,iphi,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(1.2);
      h_mMass[KEY_QA]->SetLineColor(kGray+2);
      h_mMass[KEY_QA]->SetTitle(Form("%d#pi/6<#beta<%d#pi/6, %1.1f/4.5<cos(#theta)<%1.1f/4.5",iphi,iphi+1,float(ithetah)-4.5,float(ithetah)-3.5));
      h_mMass[KEY_QA]->DrawCopy("pE");
      // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

      string pT_range = Form("[%d#pi/6,%d#pi/6)",i_phipsi,i_phipsi+1);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
    c_pT_rebin->Update();
    c_pT_rebin->Print(outputname.c_str());
  }
  output_stop = Form("%s]",outputname.c_str());
  c_pT_rebin->Print(output_stop.c_str()); // close pdf file
#endif


  // write background subtracted histograms to output file
  string outputfile = Form("../output/AuAu%s/%s/%s2DPhiPsi_InvMassSubBg_%s_phipsi_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str(),i_phipsi);
  //string outputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/%s2DPhiPsi_InvMassSubBg_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_pt = 0; i_pt < vmsa::pt_rebin_y_2D; i_pt++)
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  for(int i_phi = vmsa::CTS_start; i_phi < 12/*2*vmsa::CTS_stop*/; i_phi++) // cos(theta*) loop
	  {
	    for(int i_thetah = vmsa::CTS_start; i_thetah < 9/*2*vmsa::CTS_stop*/; i_thetah++) // cos(theta*) loop
	    {
	      for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
	      {
	        string KEY = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	        h_mMass[KEY]->Write();
	      }
	    }
          }
	}
      }
    }
  }
  File_OutPut->Close();
}
