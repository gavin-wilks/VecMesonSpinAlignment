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

void compRhoEff_F(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "AccRes")
{
  string inputfileF0 = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_0.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  string inputfileF1 = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_1.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());

  TFile *File_InPutF0 = TFile::Open(inputfileF0.c_str());
  TFile *File_InPutF1 = TFile::Open(inputfileF1.c_str());

  TString StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TGraphAsymmErrors *g_StatF0   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF0->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_SysF0    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF0->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_StatF1   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF1->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_SysF1    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF1->Get( SysErrorRho.Data()))->Clone();

  TH1F *h_frame = (TH1F*)File_InPutF1->Get("h_frame");

  TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

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

  besi19_sys->SetLineColor(kBlack);
  besi19_sys->SetMarkerStyle(20);
  besi19_sys->SetMarkerSize(1.3);
  besi19_sys->SetMarkerColor(kBlack);
  besi19_sys->SetLineColor(kBlack);
  besi19_sys->Draw("pE [] same");

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

  string leg_besi = "#phi BES-I";
  Draw_TGAE_Point_new_Symbol(0.5,0.49,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  plotTopLegend((char*)leg_besi.c_str(),0.6,0.487,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.47,0.0,0.0,0.0,0.0,20,kBlue,1.3); 
  string leg_count = Form("#phi F from BESII");
  plotTopLegend((char*)leg_count.c_str(),0.6,0.467,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  string leg_count3D = Form("#phi F from BESI");
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

  g_StatF0->SetPoint(0,0,-1.0);
  g_StatF1->SetPoint(0,0,-1.0);
  g_SysF0->SetPoint(0,0,-1.0);
  g_SysF1->SetPoint(0,0,-1.0);
  g_StatF0->SetPoint(1,0,-1.0);
  g_StatF1->SetPoint(1,0,-1.0);
  g_SysF0->SetPoint(1,0,-1.0);
  g_SysF1->SetPoint(1,0,-1.0);
  
  g_StatF0->SetMarkerStyle(20);
  g_StatF0->SetMarkerColor(kBlue);
  g_StatF0->SetLineColor(kBlue);
  g_StatF0->SetMarkerSize(1.3);

  g_SysF0->SetMarkerStyle(20);
  g_SysF0->SetMarkerSize(1.3);
  g_SysF0->SetMarkerColor(kBlue);
  g_SysF0->SetLineColor(kGray+2);
  g_StatF0->Draw("pE Z same");

  plotSysErrorsBox(g_SysF0,kBlue,energy);

  g_StatF1->SetMarkerStyle(20);
  g_StatF1->SetMarkerColor(kOrange+7);
  g_StatF1->SetLineColor(kOrange+7);
  g_StatF1->SetMarkerSize(1.3);

  g_SysF1->SetMarkerStyle(20);
  g_SysF1->SetMarkerSize(1.3);
  g_SysF1->SetMarkerColor(kOrange+7);
  g_SysF1->SetLineColor(kOrange+7);
  g_StatF1->Draw("pE Z same");

  plotSysErrorsBox(g_SysF1,kOrange+7,energy);
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
