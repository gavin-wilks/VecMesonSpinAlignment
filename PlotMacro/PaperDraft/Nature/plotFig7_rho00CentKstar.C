#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../../Utility/draw.h"
// #include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotFig7_rho00CentKstar(int energy = 6)
{
  gStyle->SetOptDate(0);
  const int style_Kstr = 20;
  const int color_Kstr = kAzure+2;

  const float size_marker = 1.4;
  
  int const NumBeamEnergy = 7;
  std::string const mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV_Run18","39GeV","54GeV","200GeV"};
  std::string const mLegEnergy[NumBeamEnergy] = {"7 GeV","11 GeV","19 GeV","27 GeV (Run18)","39 GeV","54.4 GeV","200 GeV"};

  string inputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  string GraName_2nd_stat = Form("g_rhoCent_%s_2nd_stat",mBeamEnergy[energy].c_str());
  TGraphAsymmErrors *g_rhoCent_2nd_stat = (TGraphAsymmErrors*)File_InPut->Get(GraName_2nd_stat.c_str());

  string GraName_2nd_sys = Form("g_rhoCent_%s_2nd_sys",mBeamEnergy[energy].c_str());
  TGraphAsymmErrors *g_rhoCent_2nd_sys = (TGraphAsymmErrors*)File_InPut->Get(GraName_2nd_sys.c_str());

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,0,100);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  // c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,80.0);
  h_frame->GetXaxis()->SetNdivisions(510,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("Centrality (%)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.18,0.38);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  // h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,80.0,1.0/3.0,1.0/3.0,1,3,2);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoCent_2nd_stat,style_Kstr,color_Kstr,size_marker);
  plotSysErrors(g_rhoCent_2nd_sys,color_Kstr);

  Draw_TGAE_Point_new_Symbol(10,0.35,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker);
  plotTopLegend((char*)"K^{*0} (2^{nd}-order EP)",13,0.3475,0.04,1,0.0,42,0);

  string leg_energy = Form("AuAu %s (|y| < 0.5)", mLegEnergy[energy].c_str());
  plotTopLegend((char*)leg_energy.c_str(),0.38,0.25,0.04,1,0.0,42,1);
  plotTopLegend((char*)"1.0 < p_{T}< 1.5 GeV/c",0.38,0.20,0.04,1,0.0,42,1);

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig7_rhoCentKstar_%s.eps",mBeamEnergy[energy].c_str());
  c_rho00->SaveAs(FigName.c_str());
  FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig7_rhoCentKstar_%s.png",mBeamEnergy[energy].c_str());
  c_rho00->SaveAs(FigName.c_str());
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);

    PlotLine(energy-1,energy+1,rho+err,rho+err,plot_color,2,1);
    PlotLine(energy-1,energy-1,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy+1,energy+1,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy-1,energy+1,rho-err,rho-err,plot_color,2,1);
    PlotLine(energy-1,energy-1,rho-err+0.001,rho-err,plot_color,2,1);
    PlotLine(energy+1,energy+1,rho-err+0.001,rho-err,plot_color,2,1);
  }
}
