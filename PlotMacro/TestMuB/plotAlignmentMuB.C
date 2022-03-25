#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>

#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void plotAlignmentMuB()
{
  TFile *File_InPutLambda = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/TestMuB/Lambda/PHMuB.root");
  TGraphAsymmErrors *g_PHMuB_stat = (TGraphAsymmErrors*)File_InPutLambda->Get("g_PHMuB_stat")->Clone();
  TGraphAsymmErrors *g_PHMuB_sys = (TGraphAsymmErrors*)File_InPutLambda->Get("g_PHMuB_sys")->Clone();

  TFile *File_InPutVecMeson = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/TestMuB/VecMeson/rho00MuB.root");
  TGraphAsymmErrors *g_rhoMuBPhi_stat = (TGraphAsymmErrors*)File_InPutVecMeson->Get("g_rhoMuBPhi_stat")->Clone();
  TGraphAsymmErrors *g_rhoMuBPhi_sys = (TGraphAsymmErrors*)File_InPutVecMeson->Get("g_rhoMuBPhi_sys")->Clone();
  TGraphAsymmErrors *g_rhoMuBKstar_stat = (TGraphAsymmErrors*)File_InPutVecMeson->Get("g_rhoMuBKstar_stat")->Clone();
  TGraphAsymmErrors *g_rhoMuBKstar_sys = (TGraphAsymmErrors*)File_InPutVecMeson->Get("g_rhoMuBKstar_sys")->Clone();

  TCanvas *c_alignment = new TCanvas("c_alignment","c_alignment",10,10,800,800);
  c_alignment->Divide(1,2,0,0);
  for(int iPad = 0; iPad < 2; ++iPad)
  {
    c_alignment->cd(iPad+1);
    c_alignment->cd(iPad+1)->SetLeftMargin(0.15);
    c_alignment->cd(iPad+1)->SetTicks(1,1);
    c_alignment->cd(iPad+1)->SetGrid(0,0);
  }
  c_alignment->cd(1)->SetTopMargin(0.1);
  c_alignment->cd(2)->SetBottomMargin(0.20);

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,850.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetTitle("Baryon Chemical Potential #mu_{B} (MeV)");
  h_frame->GetXaxis()->SetTitleSize(0.08);
  h_frame->GetXaxis()->SetTitleOffset(0.9);
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.06);

  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitleSize(0.08);
  h_frame->GetYaxis()->SetTitleOffset(0.8);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.06);
  h_frame->DrawCopy("pE");

  const int style_Lambda = 21;
  const int color_Lambda = kGray+2;

  const int style_phi = 29;
  const int color_phi = kRed-4;
  // const int style_phi_ALICE = 30;
  // const int color_phi_ALICE = kGray+1;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;
  // const int style_Kstr_ALICE = 24;
  // const int color_Kstr_ALICE = kGray+1;
  const float size_marker = 1.4;

  c_alignment->cd(1);
  h_frame->GetYaxis()->SetRangeUser(-0.5,6.1);
  h_frame->GetYaxis()->SetTitle("Spin Polarization P_{H} (%)");
  h_frame->DrawCopy("PE");
  PlotLine(0.0,850.0,0.0,0.0,1,2,9);
  g_PHMuB_stat->SetMarkerStyle(style_Lambda);
  g_PHMuB_stat->SetMarkerColor(color_Lambda);
  g_PHMuB_stat->SetMarkerSize(1.2);
  g_PHMuB_stat->Draw("pE same");
  plotSysErrorsBox(g_PHMuB_sys,color_Lambda);

  TLegend *leg_Lambda = new TLegend(0.20,0.75,0.35,0.90);
  leg_Lambda->SetBorderSize(0);
  leg_Lambda->SetFillColor(10);
  leg_Lambda->SetFillStyle(0);
  leg_Lambda->AddEntry(g_PHMuB_stat,"#Lambda","P");
  leg_Lambda->Draw("same");

  c_alignment->cd(2);
  h_frame->GetYaxis()->SetRangeUser(0.27,0.43);
  h_frame->GetYaxis()->SetTitle("Spin Alignment #rho_{00}");
  h_frame->DrawCopy("PE");
  PlotLine(0.0,850.0,1.0/3.0,1.0/3.0,1,2,9);

  g_rhoMuBKstar_stat->SetMarkerStyle(style_Kstr);
  g_rhoMuBKstar_stat->SetMarkerColor(color_Kstr);
  g_rhoMuBKstar_stat->SetMarkerSize(1.2);
  g_rhoMuBKstar_stat->Draw("pE same");
  plotSysErrorsBox(g_rhoMuBKstar_sys,color_Kstr);

  g_rhoMuBPhi_stat->SetMarkerStyle(style_phi);
  g_rhoMuBPhi_stat->SetMarkerColor(color_phi);
  g_rhoMuBPhi_stat->SetMarkerSize(1.8);
  g_rhoMuBPhi_stat->Draw("pE same");
  plotSysErrorsBox(g_rhoMuBPhi_sys,color_phi);

  TLegend *leg_vec = new TLegend(0.20,0.80,0.35,0.95);
  leg_vec->SetBorderSize(0);
  leg_vec->SetFillColor(10);
  leg_vec->SetFillStyle(0);
  leg_vec->AddEntry(g_rhoMuBPhi_stat,"#phi","P");
  leg_vec->AddEntry(g_rhoMuBKstar_stat,"K^{*0}","P");
  leg_vec->Draw("same");

  c_alignment->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/TestMuB/alignmentMuB.eps");
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);

    PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,plot_color,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,plot_color,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,plot_color,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,plot_color,2,1);
  }
}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double errHigh = g_rho->GetErrorYhigh(i_energy);
    double errLow  = g_rho->GetErrorYlow(i_energy);
    // cout << "errLow = " << errLow << ", errHigh = " << errHigh << endl;

    // bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(energy-7.5,rho-errLow,energy+7.5,rho+errHigh);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    if(errHigh > 0 || errLow > 0) bSys[i_energy]->Draw("l Same");
  }
}
