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

void compRhoEff_3DRandom(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Eff")
{
  string inputfile = Form("../output/AuAu%s/%s/Rho_AccResSysErrors_F_0.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile3D = Form("../output/AuAu%s/%s/3DRandom/Rho_%sSysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TFile *File_InPut3D = TFile::Open(inputfile3D.c_str());

  TString StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TGraphAsymmErrors *g_Stat   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPut->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_Sys    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPut->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_Stat3D = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut3D->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_Sys3D  = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut3D->Get( SysErrorRho.Data()))->Clone();

  TH1F *h_frame = (TH1F*)File_InPut->Get("h_frame");

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

  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  //TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

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

  Draw_TGAE_Point_new_Symbol(0.5,0.47,0.0,0.0,0.0,0.0,20,kBlue,1.3); 
  string leg_count = Form("#phi 2^{nd} order EP");
  plotTopLegend((char*)leg_count.c_str(),0.6,0.467,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  string leg_count3D = Form("#phi 3D random EP");
  plotTopLegend((char*)leg_count3D.c_str(),0.6,0.447,0.03,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  //string leg_countA = "Eff+Acc+Res Corrected #phi BES-II";
  //plotTopLegend((char*)leg_countA.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  PlotLine(0.25,0.75,0.262,0.262,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),0.85,0.26,0.03,1,0.0,42,0);


  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.26,0.04,1,0.0,42,0);
  plotTopLegend((char*)"20%-60%",3.1,0.24,0.04,1,0.0,42,0);

  g_Stat->SetPoint(0,0,-1.0);
  g_Stat3D->SetPoint(0,0,-1.0);
  g_Sys->SetPoint(0,0,-1.0);
  g_Sys3D->SetPoint(0,0,-1.0);
  g_Stat->SetPoint(1,0,-1.0);
  g_Stat3D->SetPoint(1,0,-1.0);
  g_Sys->SetPoint(1,0,-1.0);
  g_Sys3D->SetPoint(1,0,-1.0);
  
  g_Stat->SetMarkerStyle(20);
  g_Stat->SetMarkerColor(kBlue);
  g_Stat->SetLineColor(kBlue);
  g_Stat->SetMarkerSize(1.3);

  g_Sys->SetMarkerStyle(20);
  g_Sys->SetMarkerSize(1.3);
  g_Sys->SetMarkerColor(kBlue);
  g_Sys->SetLineColor(kGray+2);
  g_Stat->Draw("pE Z same");

  plotSysErrorsBox(g_Sys,kBlue,energy);

  g_Stat3D->SetMarkerStyle(20);
  g_Stat3D->SetMarkerColor(kOrange+7);
  g_Stat3D->SetLineColor(kOrange+7);
  g_Stat3D->SetMarkerSize(1.3);

  g_Sys3D->SetMarkerStyle(20);
  g_Sys3D->SetMarkerSize(1.3);
  g_Sys3D->SetMarkerColor(kOrange+7);
  g_Sys3D->SetLineColor(kOrange+7);
  g_Stat3D->Draw("pE Z same");

  plotSysErrorsBox(g_Sys3D,kOrange+7,energy);
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
