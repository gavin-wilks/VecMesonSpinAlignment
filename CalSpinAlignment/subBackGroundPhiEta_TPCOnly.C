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

void subBackGroundPhiEta_TPCOnly(int energy = 4, int i_eta = 0, int pid = 0, int year = 0, string date = "20241003", bool random3D = false, int order = 2, std::string etamode = "eta1p5_eta1p5")
{
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  std::string EP[2] = {"1st","2nd"};

  string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_TPCOnly.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(order == 1) InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(order == 1 && random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_%s_%s_%s_FirstOrder_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //if(energy == 3) InPutFile_SE = Form("../data/Yields_Phi_SE_14GeV_%s.root",etamode.c_str());
  //string InPutFile_SE = Form("../data/Yields_Phi_SE_19GeV_20220910_eta1.root");
  //if(random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_19GeV_%s_3DRandom.root",date.c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_TPCOnly.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(order == 1) InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(order == 1 && random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_%s_%s_%s_FirstOrder_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //if(energy == 3) InPutFile_ME = Form("../data/Yields_Phi_ME_14GeV_%s.root",etamode.c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_19GeV_20220910_eta1.root");
  //if(random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_19GeV_%s_3DRandom.root",date.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
  {
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
    {
      //for(Int_t i_eta = 0; i_eta < vmsa::eta_total; i_eta++) // centrality bin
      //{
        for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
        {
          for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
          {
            for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
            {
              if( i_dca != 0 && i_sig != 0 ) continue;
              string KEY_InPutSE = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_SE",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
              h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 

              string KEY_InPutME = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_ME",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
              h_mInPut_ME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 

              for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
              {
                string KEY_SE = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());
                string KEY_ME = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_ME.c_str());
                string KEY_SM = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
                if(i_norm < 2)
                {
          	  int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_norm]);
          	  int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_norm]);

          	  float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
          	  float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);

                  //cout << "Inte_SE = " << Inte_SE << "     Inte_ME = " << Inte_ME << endl;
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
                    float iSE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
                    float iME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
          	    Inte_SE += iSE;
          	    Inte_ME += iME;
                    //cout << "iSE = " << iSE << "      iME = " << iME << endl;
          	  }

                  //cout << "Inte_SE = " << Inte_SE << "     Inte_ME = " << Inte_ME << endl;
          	  if(Inte_ME != 0.0) h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
          	  h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
                }
                //cout << KEY_SM << endl;
              }
            }
          }
        }
     // }
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
  string KEY_SE_QA = Form("eta_%d_pt_%d_Centrality_1_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_eta,7,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("eta_%d_pt_%d_Centrality_1_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_eta,7,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("eta_%d_pt_%d_Centrality_1_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_eta,7,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
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


  // QA Plots for pT bins
//  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
//  c_pT->Divide(5,5);
//  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
//  {
//    c_pT->cd(i_pt+1);
//    c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
//    c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
//    c_pT->cd(i_pt+1)->SetTicks(1,1);
//    c_pT->cd(i_pt+1)->SetGrid(0,0);
//    string KEY_SE_QA = Form("eta_%d_pt_%d_Centrality_1_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_eta,i_pt,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//    h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
//    h_mMass_SE[KEY_SE_QA]->DrawCopy();
//
//    string KEY_ME_QA = Form("eta_%d_pt_%d_Centrality_1_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_eta,i_pt,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
//    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
//    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
//    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");
//
//    string KEY_SM_QA = Form("eta_%d_pt_%d_Centrality_1_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_eta,i_pt,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//    h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
//    h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
//    h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
//    h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");
//
//    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
//    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
//
//    if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
//    {
//      PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
//      PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
//    }
//    if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
//    {
//      PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
//      PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
//    }
//  }
  // QA Plots for pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
  c_pT->Divide(5,5);

  string outputname = Form("figures/%s/%s/rapiditystudy/TPCOnly_InvMassDistributions_%s_Order%d_y%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_eta);
  string output_start = Form("%s[",outputname.c_str());
  
  c_pT->Print(output_start.c_str());
  for(int itheta = 0; itheta < 7; itheta++)
  {
    for(int icent = 0; icent < 9; icent++)
    {
      for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
      {
        c_pT->cd(i_pt+1);
        c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
        c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
        c_pT->cd(i_pt+1)->SetTicks(1,1);
        c_pT->cd(i_pt+1)->SetGrid(0,0);
        string KEY_SE_QA = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_eta,i_pt,icent,itheta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        //cout << KEY_SE_QA << endl;
        h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
        h_mMass_SE[KEY_SE_QA]->SetTitle(Form("%.1f<y<%.1f, Cent %d, %d/7<|cos(#theta*)|<%d/7",vmsa::eta_bin[i_eta],vmsa::eta_bin[i_eta+1],icent,itheta,itheta+1));
        h_mMass_SE[KEY_SE_QA]->DrawCopy();

        string KEY_ME_QA = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_eta,i_pt,icent,itheta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        //cout << KEY_ME_QA << endl;
        h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
        h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
        h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
        h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

        string KEY_SM_QA = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_eta,i_pt,icent,itheta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        //cout << KEY_SM_QA << endl;
        h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
        h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
        h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
        h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

        string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
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
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_pT->Print(output_stop.c_str()); // close pdf file
//#endif
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  //for(int i_eta = 0; i_eta < vmsa::eta_total; i_eta++) // Centrality loop
  //{
   // for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
   // {
      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
      {
        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
          {
            if( i_dca != 0 && i_sig != 0 ) continue;
            for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
            {
              for(int pt_bin = vmsa::pt_rebin_first_y[energy]; pt_bin <= vmsa::pt_rebin_last_y[energy]; pt_bin++) // pt loop
              {
                //cout << "pt bin " << pt_bin << " = [" << vmsa::ptRawStart[vmsa::pt_rebin_start_y[energy][pt_bin]] << ", " << vmsa::ptRawStop[vmsa::pt_rebin_stop_y[energy][pt_bin]] << "]" << endl;
                for(int cent_bin = 0; cent_bin < vmsa::cent_rebin_total_2060; cent_bin++) // Centrality loop
                {
                  string KEY = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,pt_bin,cent_bin,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  for(int i_pt = vmsa::pt_rebin_start_y[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop_y[energy][pt_bin]; i_pt++)
                  {
                    //cout << "ipt = [" << vmsa::ptRawStart[i_pt] << ", " << vmsa::ptRawStop[i_pt] << "]" << endl;
                    for(int i_cent = vmsa::cent_rebin_2060[cent_bin]; i_cent < vmsa::cent_rebin_2060[cent_bin+1]; i_cent++)
                    {
                      //cout << "pT = " << vmsa::pt_bin[i_pt] << ", ";
                      string KEY_SM = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                      if(i_pt == vmsa::pt_rebin_start_y[energy][pt_bin] && i_cent == vmsa::cent_rebin_2060[cent_bin])
                      {
                        h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone(KEY.c_str());
                        //cout << "Create " << KEY << " with " << KEY_SM << endl;
                      }
                      else
                      {
                        h_mMass[KEY]->Add(h_mMass_SM[KEY_SM],1.0);
                        //cout << "Add to " << KEY << " with " << KEY_SM << endl;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    //}
  //}

#if _PlotQA_
  // QA Plots for pT rebins
  /*TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::centStart,vmsa::CTS_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(1.2);
    h_mMass[KEY_QA]->SetLineColor(kGray+2);
    h_mMass[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }*/
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);

  outputname = Form("figures/%s/%s/rapiditystudy/TPCOnly_InvMassBackgroundSubtracted_Rebinned_%s_Order%d_y%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_eta);
  output_start = Form("%s[",outputname.c_str());
  
  c_pT_rebin->Print(output_start.c_str());

  int centlow[3] = {40,10,0};
  int centhigh[3] = {80,40,10};

  for(int itheta = 0; itheta < 7; itheta++)
  {
    for(int icent = 0; icent < vmsa::cent_rebin_total_2060; icent++)
    {
      for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; i_pt++) // pt loop
      {
        c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
        c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
        c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
        c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
        c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
        string KEY_QA = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_pt,icent,itheta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        //cout << KEY_QA << endl;
        h_mMass[KEY_QA]->SetMarkerStyle(24);
        h_mMass[KEY_QA]->SetMarkerSize(1.2);
        h_mMass[KEY_QA]->SetLineColor(kGray+2);
        h_mMass[KEY_QA]->SetTitle(Form("%.1f<y<%.1f, Cent (%d-%d), %d/7<|cos(#theta*)|<%d/7",vmsa::eta_bin[i_eta],vmsa::eta_bin[i_eta+1],centlow[icent],centhigh[icent],itheta,itheta+1));
        h_mMass[KEY_QA]->DrawCopy("pE");
        // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

        string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low_y[energy][i_pt],vmsa::pt_up_y[energy][i_pt]);
        plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
        PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
      }
      c_pT_rebin->Update();
      c_pT_rebin->Print(outputname.c_str());
    }
  }
  output_stop = Form("%s]",outputname.c_str());
  c_pT_rebin->Print(output_stop.c_str()); // close pdf file
#endif

  // Poly + Breit Wignar fit to phi integrated InvMass
  /*TH1FMap h_mMass_theta;
  vecFMap ParFit_theta;
  for(int i_eta = 0; i_eta < vmsa::eta_total; i_eta++) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      { 
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          string KEY_theta = Form("eta_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
          for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
          {
            string KEY = Form("eta_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
    	    if(i_theta == vmsa::CTS_start) h_mMass_theta[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone(KEY_theta.c_str());
    	    else h_mMass_theta[KEY_theta]->Add(h_mMass[KEY],1.0);
          }

          TF1 *f_sig = new TF1("f_sig",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
          f_sig->SetParameter(0,vmsa::InvMass[pid]);
          f_sig->SetParLimits(0,vmsa::InvMass[pid]-1.5*vmsa::Width[pid],vmsa::InvMass[pid]+1.5*vmsa::Width[pid]);
          f_sig->SetParameter(1,vmsa::Width[pid]);
          f_sig->SetParLimits(1,0.004,0.070);
          f_sig->SetParameter(2,1.0);
          f_sig->SetParameter(3,-1.0);
          f_sig->SetParameter(4,1.0);
          f_sig->SetParameter(2,h_mMass_theta[KEY_theta]->GetMaximum()/f_sig->GetMaximum());
          f_sig->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
          h_mMass_theta[KEY_theta]->Fit(f_sig,"MNR");

          ParFit_theta[KEY_theta].clear();
          for(int i_par = 0; i_par < 5; ++i_par)
          {
    	    ParFit_theta[KEY_theta].push_back(static_cast<float>(f_sig->GetParameter(i_par)));
          }
        }
      }
    }
  }*/

#if _PlotQA_
  // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
  /*TCanvas *c_mMass_theta = new TCanvas("c_mMass_theta","c_mMass_theta",10,10,1400,1400);
  c_mMass_theta->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::centStart,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_theta[KEY_theta_QA]->SetMarkerColor(1);
    h_mMass_theta[KEY_theta_QA]->SetMarkerStyle(24);
    h_mMass_theta[KEY_theta_QA]->SetMarkerSize(0.8);
    h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");

    TF1 *f_sig = new TF1("f_sig",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
    for(int i_par = 0; i_par < 5; ++i_par)
    {
      f_sig->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par]);
    }
    f_sig->SetLineColor(2);
    f_sig->SetLineStyle(1);
    f_sig->SetLineWidth(2);
    f_sig->DrawCopy("l same");

    TF1 *f_bg = new TF1("f_bg",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
    for(int i_par = 0; i_par < 2;++i_par)
    {
      f_bg->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par+3]);
    }
    f_bg->SetLineColor(4);
    f_bg->SetLineStyle(2);
    f_bg->SetLineWidth(4);
    f_bg->DrawCopy("l same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }*/
#endif

    // Poly+bw fits for phi differential InvMass
    /*vecFMap ParFit;
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
            if( i_dca != 0 && i_sig != 0 ) continue;
	    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	    {
	      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	      {
		string KEY_theta = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
		TF1 *f_sig = new TF1("f_sig",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
		f_sig->FixParameter(0,ParFit_theta[KEY_theta][0]);
		f_sig->FixParameter(1,ParFit_theta[KEY_theta][1]);
		f_sig->SetParameter(2,ParFit_theta[KEY_theta][2]/7.0);
		f_sig->SetParameter(3,ParFit_theta[KEY_theta][3]/7.0);
		f_sig->SetParameter(4,ParFit_theta[KEY_theta][4]);
		f_sig->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);

		cout << "i_pt = " << i_pt << ", i_cent = " << i_cent << ", i_theta = " << i_theta << ", i_dca = " << i_dca << ", i_sig = " << i_sig << ", i_norm = " << i_norm << endl;
		string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
		TFitResultPtr result = h_mMass[KEY]->Fit(f_sig,"MNRIS");
		ParFit[KEY].clear();
		for(int i_par = 0; i_par < 5; ++i_par)
		{
		  ParFit[KEY].push_back(static_cast<float>(f_sig->GetParameter(i_par)));
		}

		TF1 *f_bg = new TF1("f_bg",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
		for(int i_par = 0; i_par < 2;++i_par)
		{
		  f_bg->SetParameter(i_par,ParFit[KEY][i_par+3]);
		  f_bg->SetParError(i_par,f_sig->GetParError(i_par+3));
		}
                
	        //f_bg->SetParError(0,f_bw->GetParError(3));
	        //f_bg->SetParError(1,f_bw->GetParError(4));

                double params[2] = {result->GetParams()[3],result->GetParams()[4]};
                TMatrixDSym covArr(2);
                covArr(0,0) = result->GetCovarianceMatrix()(3,3);
                covArr(0,1) = result->GetCovarianceMatrix()(3,4);
                covArr(1,0) = result->GetCovarianceMatrix()(4,3);
                covArr(1,1) = result->GetCovarianceMatrix()(4,4);

                TH1F *hTotal = (TH1F*)h_mMass[KEY]->Clone();
	        h_mMass[KEY]->Add(f_bg,-1.0); // subtract linear background for phi differential InvMass
                for(int bin = h_mMass[KEY]->FindBin(vmsa::BW_Start[pid]); bin <= h_mMass[KEY]->FindBin(vmsa::BW_Stop[pid]); bin++)
                {
                  float binWidth = h_mMass[KEY]->GetBinWidth(bin);  
                  float binLowEdge = hTotal->GetBinLowEdge(bin);
                  float funcError = f_bg->IntegralError(binLowEdge,binLowEdge+binWidth,params,covArr.GetMatrixArray())/binWidth;
                  float totError = hTotal->GetBinError(bin);
                  if(i_cent == 9)
                  { 
                    cout << "funcError  = " << funcError << endl;
                    cout << "totError   = " << totError << endl;
                    cout << "finalError = " << TMath::Sqrt(totError*totError+funcError*funcError) << endl;
                  }
                  h_mMass[KEY]->SetBinError(bin,TMath::Sqrt(totError*totError+funcError*funcError)); 
                } 
                delete hTotal;

		//h_mMass[KEY]->Add(f_bg,-1.0);
	      }
	    }
	  }
	}
      }
    }*/

#if _PlotQA_
    // QA plots for phi differential InvMass after linear background subtraction
   /* TCanvas *c_mMass_phi_diff = new TCanvas("c_mMass_phi_diff","c_mMass_phi_diff",10,10,1400,1400);
    c_mMass_phi_diff->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_QA = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,vmsa::centStart,vmsa::CTS_start,vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("PE");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }*/
#endif

  // write background subtracted histograms to output file
  string outputfile = Form("../output/AuAu%s/%s/TPCOnly_InvMassSubBgEta%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),i_eta,etamode.c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/InvMassSubBgEta%d_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),i_eta,etamode.c_str());
  if(order == 1 && random3D) outputfile = Form("../output/AuAu%s/%s/3Drandom/InvMassSubBgEta%d_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),i_eta,etamode.c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  //for(int i_eta = 0; i_eta < vmsa::eta_total; i_eta++) // Centrality loop
  //{
    for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; i_pt++) // pt loop
    {
      //for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
      for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2060; i_cent++) // Centrality loop
      {
        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
          {
            if( i_dca != 0 && i_sig != 0 ) continue;
            for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
            {
              for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
              {
                string KEY = Form("eta_%d_pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                h_mMass[KEY]->Write();
              }
            }
          }
        }
      }
    }
 // }
  File_OutPut->Close();
}
