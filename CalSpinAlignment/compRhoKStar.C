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

void compRhoKStar(Int_t energy = 4, Int_t pid = 2, Int_t i_cent = 9)
{
  string inputfileRaw = Form("../output/AuAu%s/%s/RawRhoPt/RawRhoPtSys_TPCorToF_20220411.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfileEff = Form("../output/AuAu%s/%s/RawRhoPt/RawRhoPtSys_TPCandToF_20220523.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TFile *File_InPutRaw = TFile::Open(inputfileRaw.c_str());
  TFile *File_InPutEff = TFile::Open(inputfileEff.c_str());

  TString StatErrorRho = "Rho00_KStar_Centrality_9_BW";

  TGraphAsymmErrors *g_RawStat = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutRaw->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_EffStat = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPutEff->Get(StatErrorRho.Data()))->Clone();
  cout << "Loaded in graphs" << endl;
  
  
  //TH1F *h_frame = (TH1F*)File_InPutRaw->Get("h_frame");

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

  h_frame->GetYaxis()->SetRangeUser(0.2,0.55);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
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

  Draw_TGAE_Point_new_Symbol(2.5,0.47,0.0,0.0,0.0,0.0,20,kBlue,1.3); 
  string leg_count = "TPC || ToF K*^{0} BES-II";
  plotTopLegend((char*)leg_count.c_str(),2.6,0.467,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(2.5,0.45,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  string leg_countE = "TPC && ToF K*^{0} BES-II";
  plotTopLegend((char*)leg_countE.c_str(),2.6,0.447,0.03,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  //string leg_countA = "Eff+Acc+Res Corrected #phi BES-II";
  //plotTopLegend((char*)leg_countA.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  PlotLine(2.25,2.5,0.428,0.428,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),2.6,0.427,0.03,1,0.0,42,0);


  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),0.5,0.26,0.04,1,0.0,42,0);
  plotTopLegend((char*)"20%-60%",1.1,0.24,0.04,1,0.0,42,0);
  
  g_RawStat->SetMarkerStyle(20);
  g_RawStat->SetMarkerColor(kBlue);
  g_RawStat->SetLineColor(kBlue);
  g_RawStat->SetMarkerSize(1.3);

  //g_RawSys->SetMarkerStyle(20);
  //g_RawSys->SetMarkerSize(1.3);
  //g_RawSys->SetMarkerColor(kGray+2);
  //g_RawSys->SetLineColor(kGray+2);
  g_RawStat->Draw("pE Z same");

  //plotSysErrorsBox(g_RawSys,kGray+2,energy);

  g_EffStat->SetMarkerStyle(20);
  g_EffStat->SetMarkerColor(kOrange+7);
  g_EffStat->SetLineColor(kOrange+7);
  g_EffStat->SetMarkerSize(1.3);

  //g_EffSys->SetMarkerStyle(20);
  //g_EffSys->SetMarkerSize(1.3);
  //g_EffSys->SetMarkerColor(kBlue);
  //g_EffSys->SetLineColor(kRed);
  g_EffStat->Draw("pE Z same");

  //plotSysErrorsBox(g_EffSys,kBlue,energy);

  //g_AccStat->SetMarkerStyle(20);
  //g_AccStat->SetMarkerColor(kOrange+7);
  //g_AccStat->SetLineColor(kOrange+7);
  //g_AccStat->SetMarkerSize(1.3);

  //g_AccSys->SetMarkerStyle(20);
  //g_AccSys->SetMarkerSize(1.3);
  //g_AccSys->SetMarkerColor(kOrange+7);
  //g_AccSys->SetLineColor(kOrange+7);
  //g_AccStat->Draw("pE Z same");

  //plotSysErrorsBox(g_AccSys,kOrange+7,energy);
  for (int ipt = 0; ipt < 6; ipt++)
  {
    double pt;
    double b1val, b1stat, b1sys;
    double b2val, b2stat, b2sys;
    g_RawStat->GetPoint(ipt,pt,b1val);
    b1stat = g_RawStat->GetErrorYhigh(ipt);
    g_EffStat->GetPoint(ipt,pt,b2val);
    b2stat = g_EffStat->GetErrorYhigh(ipt);

    double sigma = TMath::Abs((b2val-b1val)/TMath::Sqrt(b1stat*b1stat + b2stat*b2stat));
    std::cout << std::fixed << std::setprecision(4);
    cout << "pt = "<<  pt << ":    ToF || TPC = " << b1val << " +/- " << b1stat << " (stat)    ToF && TPC = " << b2val << " +/- " << b2stat << " (stat)    sigma = " << std::setprecision(2) << sigma << endl; 
  }
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
