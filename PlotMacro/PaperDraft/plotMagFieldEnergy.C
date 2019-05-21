#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TMath.h"

using namespace std;

float ErrorAdd(float x, float y)
{
  return sqrt(x*x+y*y);
}

float ErrTimes(float x, float y, float dx, float dy)
{
  return x*y*ErrorAdd(dx/x,dy/y);
}

float ErrDiv(float x, float y, float dx, float dy)
{
  return x/y*ErrorAdd(dx/x,dy/y);
}

void plotMagFieldEnergy()
{
  bool isPlotMean = true;

  gStyle->SetOptDate(0);
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/magFiled_energy_SysErrors.root");
  TGraphAsymmErrors *g_eMagOverT_1st     = (TGraphAsymmErrors*)File_Input->Get("g_eMagOverT_1st");
  TGraphAsymmErrors *g_eMag_1st          = (TGraphAsymmErrors*)File_Input->Get("g_eMag_1st");
  TGraphAsymmErrors *g_eMagTesla_1st     = (TGraphAsymmErrors*)File_Input->Get("g_eMagTesla_1st");
  TGraphAsymmErrors *g_eMeanMagOverT_1st = (TGraphAsymmErrors*)File_Input->Get("g_eMeanMagOverT_1st");
  TGraphAsymmErrors *g_eMeanMag_1st      = (TGraphAsymmErrors*)File_Input->Get("g_eMeanMag_1st");
  TGraphAsymmErrors *g_eMeanMagTesla_1st = (TGraphAsymmErrors*)File_Input->Get("g_eMeanMagTesla_1st");

  TGraphAsymmErrors *g_eMagOverT_2nd     = (TGraphAsymmErrors*)File_Input->Get("g_eMagOverT_2nd");
  TGraphAsymmErrors *g_eMag_2nd          = (TGraphAsymmErrors*)File_Input->Get("g_eMag_2nd");
  TGraphAsymmErrors *g_eMagTesla_2nd     = (TGraphAsymmErrors*)File_Input->Get("g_eMagTesla_2nd");
  TGraphAsymmErrors *g_eMeanMagOverT_2nd = (TGraphAsymmErrors*)File_Input->Get("g_eMeanMagOverT_2nd");
  TGraphAsymmErrors *g_eMeanMag_2nd      = (TGraphAsymmErrors*)File_Input->Get("g_eMeanMag_2nd");
  TGraphAsymmErrors *g_eMeanMagTesla_2nd = (TGraphAsymmErrors*)File_Input->Get("g_eMeanMagTesla_2nd");


  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-100.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(9.0,240.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(-10.0,300.0);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("e#bar{B}/T (MeV)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();

  //-------------------------------plot of eB/T [MeV]-------------------------------
  TCanvas *c_eMagOverT = new TCanvas("c_eMagOverT","c_eMagOverT",10,10,800,800);
  c_eMagOverT->cd();
  c_eMagOverT->cd()->SetLeftMargin(0.15);
  c_eMagOverT->cd()->SetBottomMargin(0.15);
  c_eMagOverT->cd()->SetTicks(1,1);
  c_eMagOverT->cd()->SetGrid(0,0);
  c_eMagOverT->cd()->SetLogx();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,0.0,0.0,1,3,2);

  plotTopLegend((char*)"Au+Au (20-60\% & |#eta| < 1)",0.45,0.85,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",0.4,0.80,0.04,1,0.0,42,1);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMagOverT_1st,20,kAzure+2,1.4);
  Draw_TGAE_Point_new_Symbol(70,230,0.0,0.0,0.0,0.0,20,kAzure+2,1.4);
  plotTopLegend((char*)"1^{st}-order EP",80,225,0.04,1,0.0,42,0);
  if(isPlotMean) Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMeanMagOverT_1st,24,kAzure+2,1.4);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMagOverT_2nd,29,kRed,1.8);
  Draw_TGAE_Point_new_Symbol(70,210,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"2^{nd}-order EP",80,205,0.04,1,0.0,42,0);
  if(isPlotMean) Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMeanMagOverT_2nd,30,kRed,1.8);

  plotTopLegend((char*)"#frac{eMag}{T} = 3m_{s}#sqrt{9(#rho_{00} - #frac{1}{3})}",30,40,0.04,1,0.0,42,0);

  c_eMagOverT->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/c_eMagOverT_energy.eps");
  //-------------------------------plot of eB/T [MeV]-------------------------------

  //-------------------------------plot of eB [(MeV)^2]-------------------------------
  TCanvas *c_eMag = new TCanvas("c_eMag","c_eMag",10,10,800,800);
  c_eMag->cd();
  c_eMag->cd()->SetLeftMargin(0.15);
  c_eMag->cd()->SetBottomMargin(0.15);
  c_eMag->cd()->SetTicks(1,1);
  c_eMag->cd()->SetGrid(0,0);
  c_eMag->cd()->SetLogx();
  c_eMag->cd()->SetLogy();
  h_frame->GetYaxis()->SetRangeUser(9.0e2,2.0e5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("e#bar{B} ((MeV)^{2})");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  // PlotLine(9.0,240.0,0.0,0.0,1,3,2);

  plotTopLegend((char*)"Au+Au (20-60\% & |#eta| < 1)",0.45,0.85,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",0.4,0.80,0.04,1,0.0,42,1);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMag_1st,20,kAzure+2,1.4);
  Draw_TGAE_Point_new_Symbol(70,6.5e4,0.0,0.0,0.0,0.0,20,kAzure+2,1.4);
  plotTopLegend((char*)"1^{st}-order EP",80,6.2e4,0.04,1,0.0,42,0);
  if(isPlotMean) Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMeanMag_1st,24,kAzure+2,1.4);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMag_2nd,29,kRed,1.8);
  Draw_TGAE_Point_new_Symbol(70,4.5e4,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"2^{nd}-order EP",80,4.2e4,0.04,1,0.0,42,0);
  if(isPlotMean) Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMeanMag_2nd,30,kRed,1.8);

  plotTopLegend((char*)"T = 120 MeV",30,3e3,0.04,1,0.0,42,0);

  c_eMag->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/c_eMag_energy.eps");
  //-------------------------------plot of eB [(MeV)^2]-------------------------------

  //-------------------------------plot of eB [e*T]-------------------------------
  TCanvas *c_eMagTesla = new TCanvas("c_eMagTesla","c_eMagTesla",10,10,800,800);
  c_eMagTesla->cd();
  c_eMagTesla->cd()->SetLeftMargin(0.15);
  c_eMagTesla->cd()->SetBottomMargin(0.15);
  c_eMagTesla->cd()->SetTicks(1,1);
  c_eMagTesla->cd()->SetGrid(0,0);
  c_eMagTesla->cd()->SetLogx();
  c_eMagTesla->cd()->SetLogy();
  h_frame->GetYaxis()->SetRangeUser(0.5e14,1e15);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("e#bar{B} (eT)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  // PlotLine(9.0,240.0,0.0,0.0,1,3,2);

  plotTopLegend((char*)"Au+Au (20-60\% & |#eta| < 1)",0.45,0.85,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",0.4,0.80,0.04,1,0.0,42,1);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMagTesla_1st,20,kAzure+2,1.4);
  Draw_TGAE_Point_new_Symbol(70,5.5e14,0.0,0.0,0.0,0.0,20,kAzure+2,1.4);
  plotTopLegend((char*)"1^{st}-order EP",80,5.2e14,0.04,1,0.0,42,0);
  if(isPlotMean) Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMeanMagTesla_1st,24,kAzure+2,1.4);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMagTesla_2nd,29,kRed,1.8);
  Draw_TGAE_Point_new_Symbol(70,4.5e14,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"2^{nd}-order EP",80,4.2e14,0.04,1,0.0,42,0);
  if(isPlotMean) Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_eMeanMagTesla_2nd,30,kRed,1.8);

  plotTopLegend((char*)"T = 120 MeV",30,3e3,0.04,1,0.0,42,0);

  c_eMagTesla->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/c_eMagTesla_energy.eps");
  //-------------------------------plot of eB [(MeV)^2]-------------------------------
}
