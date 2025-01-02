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

void subBackGround_AllOneBinFromStart(int energy = 6, int pid = 0, int year = 0, string date = "20240827", bool random3D = false, int order = 2, string etamode = "eta1_eta1")
{
  std::string EP[2] = {"","2nd"};
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //string InPutFile_SE = Form("../data/Yields_Phi_SE_19GeV_20220527.root"); //original eta < 1.0
  string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //if(energy == 3) InPutFile_SE = Form("../data/Yields_Phi_SE_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_19GeV_20220408.root"); //original eta < 1.0
  string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_NoPsi2Bin.root",vmsa::mBeamEnergy[energy].c_str(),"20240712",etamode.c_str());
  //if(energy == 3) InPutFile_ME = Form("../data/Yields_Phi_ME_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());


  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  {
    for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;

	  string KEY_InPutSE_All = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_SE",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	  string KEY_InPutME_All = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_ME",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());

          for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
          {
	    string KEY_InPutSE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_SE",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
            h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 
            //h_mInPut_SE[KEY_InPutSE]->Scale(1.89/2.443373);
            //h_mInPut_SE[KEY_InPutSE]->Rebin(2);
	    string KEY_InPutME = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_ME",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	    h_mInPut_ME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 
            //h_mInPut_ME[KEY_InPutME]->Scale(1.89/2.443373);
            //h_mInPut_ME[KEY_InPutME]->Rebin(2);
            
            if(i_pt == vmsa::pt_rebin_start[energy][0])
            {
              cout << "Start at i_pt = " << i_pt << endl;
              h_mInPut_SE[KEY_InPutSE_All] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_InPutSE_All.c_str());
              h_mInPut_ME[KEY_InPutME_All] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_InPutME_All.c_str());
            }
            if(i_pt >  vmsa::pt_rebin_start[energy][0] && i_pt <= vmsa::pt_rebin_stop[energy][0])
            {
              cout << "Add i_pt = " << i_pt << endl;
              h_mInPut_SE[KEY_InPutSE_All]->Add(h_mInPut_SE[KEY_InPutSE],1.0);
              h_mInPut_ME[KEY_InPutME_All]->Add(h_mInPut_ME[KEY_InPutME],1.0);
            }
          }
	  for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
	  {
	    string KEY_SE = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE_All]->Clone(KEY_SE.c_str());
	    string KEY_ME = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMass_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME_All]->Clone(KEY_ME.c_str());
	    string KEY_SM = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
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

  string outputfile_comp = Form("../output/AuAu%s/%s/AllRawHists_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode.c_str());

  TFile *File_OutPut_Comp = new TFile(outputfile_comp.c_str(),"RECREATE");
  File_OutPut_Comp->cd();
  for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  {
    for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;

	  string KEY_InPutSE_All = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_SE",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	  string KEY_InPutME_All = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_ME",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());

          for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
          {
	    string KEY_InPutSE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_SE",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
            if(i_cent == 9 && i_dca == 0 && i_sig == 0) 
            {
              h_mInPut_SE[KEY_InPutSE]->SetName(Form("pt_%d_Centrality_9_CosThetaStar_%d_SE",i_pt,i_theta));
              h_mInPut_SE[KEY_InPutSE]->Write();
            }        

	    string KEY_InPutME = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_ME",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
            if(i_cent == 9 && i_dca == 0 && i_sig == 0) 
            {
              h_mInPut_ME[KEY_InPutME]->SetName(Form("pt_%d_Centrality_9_CosThetaStar_%d_ME",i_pt,i_theta));
              h_mInPut_ME[KEY_InPutME]->Write();
            }        
          }
          if(i_cent == 9 && i_dca == 0 && i_sig == 0) 
          {
            h_mInPut_SE[KEY_InPutSE_All]->SetName(Form("Centrality_9_CosThetaStar_%d_SE",i_theta));
            h_mInPut_SE[KEY_InPutSE_All]->Write();
            h_mInPut_ME[KEY_InPutME_All]->SetName(Form("Centrality_9_CosThetaStar_%d_ME",i_theta));
            h_mInPut_ME[KEY_InPutME_All]->Write();
          }        
        }     
      }
    }
  }
  File_OutPut_Comp->Close();

#if _PlotQA_
  // QA Plots for SE vs. ME
    TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
    c_peak->cd();
    c_peak->cd()->SetLeftMargin(0.15);
    c_peak->cd()->SetBottomMargin(0.15);
    c_peak->cd()->SetTicks(1,1);
    c_peak->cd()->SetGrid(0,0);
    string KEY_SE_QA = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",9,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
    h_mMass_SE[KEY_SE_QA]->SetStats(0);
    h_mMass_SE[KEY_SE_QA]->SetTitle("");
    h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

    string KEY_ME_QA = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",9,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",9,vmsa::CTS_start,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
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

    c_peak->SaveAs(Form("figures/%s/%s/pTstudy/OneBigBin_CloseUpMixedEvent_%s_Order%d_Cent9.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order));

#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
	  {
	    for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	    {
	      string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",pt_bin,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              string KEY_SM = Form("Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM",i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      cout << "KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	      h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone(KEY.c_str());
	    }
	  }
	}
      }
    }
  }

  // write background subtracted histograms to output file
  string outputfile = Form("../output/AuAu%s/%s/OneBigBinFromStart_InvMassSubBg_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  //string outputfile = Form("../output/AuAu%s/%s/OneBigBinFromStart_InvMassSubBg_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  //string outputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/InvMassSubBg_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	  {
	    for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
	    {
	      string KEY = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMass[KEY]->Write();
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}
