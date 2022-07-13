#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
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

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);

void calSysError_Eff(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Eff")
{
  string inputfileHframe = Form("../output/AuAu%s/%s/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  TFile *File_InPutHframe = TFile::Open(inputfileHframe.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mRho;
  TH1DMap h_mRho;  
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
	    g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str());
	    string KEY_rho_ptbin = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	    h_mRho[KEY_rho_ptbin] = (TH1D*)File_InPut->Get(KEY_rho_ptbin.c_str());
	  }
	}
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPutHframe->Get("h_frame");

#if _PlotQA_
  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

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
	    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
	  }
	}
      }
    }
  }
#endif
 
  string KEY_Default = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
  cout << "DEFAULT: " << KEY_Default << endl;
  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();
    
  double sysErr[6][5]; // N points with 5 sources of systematics for each

  for(Int_t i_point = 0; i_point < 6; ++i_point)
  {
    double sysDca[vmsa::Dca_stop];
    double sysNSig[vmsa::nSigKaon_stop];
    double sysNorm[vmsa::Norm_stop];
    double sysSig[vmsa::Sig_stop];
    double sysMeth[vmsa::Method_stop];

    double pt_def, rho_def;
    g_mRho[KEY_Default]->GetPoint(i_point,pt_def,rho_def); 

    cout << "pt: " << pt_def << "   rho00: " << rho_def;

    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    { 
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysDca[i_dca] = rho_sys;
    }

    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,i_sig,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysNSig[i_sig] = rho_sys;
    }

    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),i_norm,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysNorm[i_norm] = rho_sys;
    }	 

    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysSig[i_sigma] = rho_sys;
    }

    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysMeth[i_method] = rho_sys;
    }

    Double_t rho_min[5] = { TMath::MinElement(vmsa::Dca_stop,sysDca),
                            TMath::MinElement(vmsa::nSigKaon_stop,sysNSig),
                            TMath::MinElement(vmsa::Norm_stop,sysNorm),
                            TMath::MinElement(vmsa::Sig_start,sysSig),
                            TMath::MinElement(vmsa::Method_stop,sysMeth)    };

    Double_t rho_max[5] = { TMath::MaxElement(vmsa::Dca_stop,sysDca),
                            TMath::MaxElement(vmsa::nSigKaon_stop,sysNSig),
                            TMath::MaxElement(vmsa::Norm_stop,sysNorm),
                            TMath::MaxElement(vmsa::Sig_start,sysSig),
                            TMath::MaxElement(vmsa::Method_stop,sysMeth)    };
  
    double SysError_rho = 0.0;
    for(int i = 0; i < 5; i++)
    {
      double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
      cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
      SysError_rho += sourcei;
    }
 
    SysError_rho = TMath::Sqrt(SysError_rho);

    Double_t pt, rho;
    g_mRho[KEY_Default]->GetPoint(i_point,pt,rho);

    //float mean_rho = total_rho/(float)counter;
    g_SysErrors->SetPoint(i_point,pt,rho);
    g_SysErrors->SetPointError(i_point,0.0,0.0,SysError_rho,SysError_rho);

    double StatError_rho = g_mRho[KEY_Default]->GetErrorYhigh(i_point);
    g_StatErrors->SetPoint(i_point,pt,rho);
    g_StatErrors->SetPointError(i_point,0.0,0.0,StatError_rho,StatError_rho);
  }
 
  cout << "Finished pT bin calculations" << endl;
 
  TFile *besi = TFile::Open("../output/AuAu19GeV/Phi/BESI/EffPhiPtSys.root");
  TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get(KEY_Default.c_str());
  //TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //cout << "All good" << endl;

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->Draw("pE");
  //g_StatErrors->SetMarkerStyle(20);
  //g_StatErrors->SetMarkerColor(kGray+2);
  //g_StatErrors->SetLineColor(2);
  //g_StatErrors->Draw("pE same");
  //g_SysErrors->SetMarkerStyle(20);
  //g_SysErrors->SetMarkerColor(kGray+2);
  //g_SysErrors->SetLineColor(kGray+2);
  //g_SysErrors->Draw("pE same");
//  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  besi19->SetLineColor(kBlack);
  besi19->SetMarkerStyle(20);
  besi19->SetMarkerSize(1.3);
  besi19->SetMarkerColor(kBlack);
  besi19->SetLineColor(kBlack);
  besi19->Draw("pE Z same");

  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->SetMarkerStyle(20);
  //besi19_sys->SetMarkerSize(1.3);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->Draw("pE [] same");

  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,kRed,1.3); 
  string leg_count;
  if(correction == "Eff") leg_count = "Efficiency Corrected #phi BES-II";
  if(correction == "Raw") leg_count = "Raw #phi BES-II";
  if(correction == "AccRes") leg_count = "Eff+Acc+Res Corrected #phi BES-II";

  plotTopLegend((char*)leg_count.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  string leg_besi = "Efficiency Corrected #phi BES-I";
  Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  PlotLine(0.25,0.75,0.262,0.262,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),0.85,0.26,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.26,0.04,1,0.0,42,0);
  plotTopLegend((char*)"20%-60%",3.1,0.24,0.04,1,0.0,42,0);

  g_StatErrors->SetMarkerStyle(20);
  g_StatErrors->SetMarkerColor(kRed);
  g_StatErrors->SetLineColor(kRed);
  g_StatErrors->SetMarkerSize(1.3);

  g_SysErrors->SetMarkerStyle(20);
  g_SysErrors->SetMarkerSize(1.3);
  g_SysErrors->SetMarkerColor(kRed);
  g_SysErrors->SetLineColor(kRed);
  //g_SysErrors->Draw("pE || same");
  g_StatErrors->Draw("pE Z same");

  plotSysErrorsBox(g_SysErrors,kRed,energy);


 // string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
 // cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
 // TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  //File_OutPut->cd();
  //h_frame->Write();
  //TString StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //g_StatErrors->SetName(StatErrorRho.Data());
  //g_StatErrors->Write();
  //TString SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //g_SysErrors->SetName(SysErrorRho.Data());
  //g_SysErrors->Write();
  //File_OutPut->Close();
}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
{
  const int nEnergy = 6;
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < 6; ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);
    
    //bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(vmsa::pt_low[beamE][i_energy],rho-err,vmsa::pt_up[beamE][i_energy],rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
}
