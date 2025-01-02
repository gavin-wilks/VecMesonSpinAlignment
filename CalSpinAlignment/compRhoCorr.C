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
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compRhoCorr(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, Int_t order = 2)
{
  std::string EP = "";
  if(order == 1) EP = "_FirstOrder";
  string inputfileRaw = Form("../output/AuAu%s/%s/Rho_RawSysErrors_F_0_eta1_eta1_PolySys%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),EP.c_str());
  string inputfileRes = Form("../output/AuAu%s/%s/Rho_RawResSysErrors_F_0_eta1_eta1_PolySys%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),EP.c_str());
  string inputfileEff = Form("../output/AuAu%s/%s/Rho_EffAccSysErrors_F_0_eta1_eta1_PolySys%s_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),EP.c_str());
  string inputfileAcc = Form("../output/AuAu%s/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys%s_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),EP.c_str());

  TFile *File_InPutRaw = TFile::Open(inputfileRaw.c_str());
  TFile *File_InPutRes = TFile::Open(inputfileRes.c_str());
  TFile *File_InPutEff = TFile::Open(inputfileEff.c_str());
  TFile *File_InPutAcc = TFile::Open(inputfileAcc.c_str());

  TString StatErrorRho = Form("g_rho00_order%d_%s_%s_StatError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString SysErrorRho = Form("g_rho00_order%d_%s_%s_SysError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TGraphAsymmErrors *g_RawStat = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutRaw->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_RawSys  = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutRaw->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_ResStat = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutRes->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_ResSys  = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutRes->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_EffStat = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutEff->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_EffSys  = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutEff->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_AccStat = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutAcc->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_AccSys  = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutAcc->Get( SysErrorRho.Data()))->Clone();

  TH1F *h_frame = (TH1F*)File_InPutRaw->Get("h_frame");

  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  //TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //cout << "All good" << endl;

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->GetYaxis()->SetRangeUser(0.3,0.42);
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

  //besi19->SetLineColor(kBlack);
  //besi19->SetMarkerStyle(20);
  //besi19->SetMarkerSize(1.3);
  //besi19->SetMarkerColor(kBlack);
  //besi19->SetLineColor(kBlack);
  //besi19->Draw("pE Z same");

  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->SetMarkerStyle(20);
  //besi19_sys->SetMarkerSize(1.3);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->Draw("pE [] same");

  TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //besi19->SetMarkerStyle(20);
  //besi19->SetMarkerSize(1.3);
  //besi19->SetMarkerColor(kBlack);
  //besi19->SetLineColor(kBlack);
  //besi19->Draw("pE Z same");

  //besi19_sys->SetMarkerStyle(20);
  //besi19_sys->SetMarkerSize(1.3);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->Draw("pE [] same");

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.49,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.487,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlue,1.3); 
  string leg_count = "Raw #phi BES-II";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.40,0.0,0.0,0.0,0.0,20,kRed,1.3); 
  string leg_countR = "Raw Resolution Corrected #phi BES-II";
  plotTopLegend((char*)leg_countR.c_str(),0.6,0.397,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.39,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  string leg_countE = "Eff x Acc Corrected #phi BES-II";
  plotTopLegend((char*)leg_countE.c_str(),0.6,0.387,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.38,0.0,0.0,0.0,0.0,20,kGray+2,1.3); 
  string leg_countA = "Fully Corrected #phi BES-II";
  plotTopLegend((char*)leg_countA.c_str(),0.6,0.377,0.03,1,0.0,42,0);

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  PlotLine(0.25,0.75,0.262,0.262,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),0.85,0.26,0.03,1,0.0,42,0);


  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.26,0.04,1,0.0,42,0);
  plotTopLegend((char*)"20%-60%",3.1,0.24,0.04,1,0.0,42,0);
  

  g_RawStat->SetMarkerStyle(20);
  g_RawStat->SetMarkerColor(kBlue);
  g_RawStat->SetLineColor(kBlue);
  g_RawStat->SetMarkerSize(1.3);

  g_RawSys->SetMarkerStyle(20);
  g_RawSys->SetMarkerSize(1.3);
  g_RawSys->SetMarkerColor(kBlue);
  g_RawSys->SetLineColor(kBlue);
  g_RawStat->Draw("pE Z same");

  plotSysErrorsBox(g_RawSys,kBlue);

  for(int i = 0; i < g_RawStat->GetN(); i++)
  {
    double pt, rho;
    g_ResStat->GetPoint(i,pt,rho);
    g_ResStat->SetPoint(i,pt+0.04,rho);
    g_ResSys->SetPoint(i,pt+0.04,rho);
  }
  g_ResStat->SetMarkerStyle(20);
  g_ResStat->SetMarkerColor(kRed);
  g_ResStat->SetLineColor(kRed);
  g_ResStat->SetMarkerSize(1.3);

  g_ResSys->SetMarkerStyle(20);
  g_ResSys->SetMarkerSize(1.3);
  g_ResSys->SetMarkerColor(kRed);
  g_ResSys->SetLineColor(kRed);
  g_ResStat->Draw("pE Z same");

  plotSysErrorsBox(g_ResSys,kRed);

  for(int i = 0; i < g_RawStat->GetN(); i++)
  {
    double pt, rho;
    g_EffStat->GetPoint(i,pt,rho);
    g_EffStat->SetPoint(i,pt+0.08,rho);
    g_EffSys->SetPoint(i,pt+0.08,rho);
  }
  g_EffStat->SetMarkerStyle(20);
  g_EffStat->SetMarkerColor(kOrange+7);
  g_EffStat->SetLineColor(kOrange+7);
  g_EffStat->SetMarkerSize(1.3);

  g_EffSys->SetMarkerStyle(20);
  g_EffSys->SetMarkerSize(1.3);
  g_EffSys->SetMarkerColor(kOrange+7);
  g_EffSys->SetLineColor(kOrange+7);
  g_EffStat->Draw("pE Z same");

  plotSysErrorsBox(g_EffSys,kOrange+7);

  for(int i = 0; i < g_RawStat->GetN(); i++)
  {
    double pt, rho;
    g_AccStat->GetPoint(i,pt,rho);
    g_AccStat->SetPoint(i,pt+0.12,rho);
    g_AccSys->SetPoint(i,pt+0.12,rho);
  }
  g_AccStat->SetMarkerStyle(20);
  g_AccStat->SetMarkerColor(kGray+2);
  g_AccStat->SetLineColor(kGray+2);
  g_AccStat->SetMarkerSize(1.3);

  g_AccSys->SetMarkerStyle(20);
  g_AccSys->SetMarkerSize(1.3);
  g_AccSys->SetMarkerColor(kGray+2);
  g_AccSys->SetLineColor(kGray+2);
  g_AccStat->Draw("pE Z same");

  plotSysErrorsBox(g_AccSys,kGray+2);
  c_rho_SysError->SaveAs(Form("./figures/Phi/%s/rho00Pt_eta1_eta1_%s_Order%d.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order));
}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
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
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.08,rho-err,pt+0.08,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}
