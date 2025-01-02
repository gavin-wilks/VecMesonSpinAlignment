#include <string>

#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>

#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotRho00Y_single()
{
  gStyle->SetOptDate(0);
  float rho00_low[6] = {0.24,0.24,0.24,0.24,0.24,0.24}; // 11.5 -- 200 GeV
  float rho00_high[6] = {0.52,0.52,0.54,0.52,0.52,0.52};
  // float rho00_low[6] = {0.20,0.20,0.20,0.24,0.24,0.29}; // 11.5 -- 200 GeV
  // float rho00_high[6] = {0.50,0.52,0.64,0.40,0.40,0.35};
  float pt_low = 0.0;
  float pt_high = 1.0;
  float font_size = 0.12;
  float leg_size = 0.10;

  string mBeanEnergy[6] = {"14.6 GeV","19.6 GeV","27 GeV","39 GeV","62.4 GeV","200 GeV"};
  int mEnergy[6] = {14,19,27,39,62,200};
  int plot_style_1st = 20;
  int plot_color_1st = kAzure+2;
  float plot_size_1st = 1.4;
  int plot_style_2nd = 29;
  int plot_color_2nd = kRed;
  float plot_size_2nd = 1.4;

  TGraphAsymmErrors *g_rho_1st_stat[6];
  TGraphAsymmErrors *g_rho_1st_sys[6];
  TGraphAsymmErrors *g_rho_2nd_stat[6];
  TGraphAsymmErrors *g_rho_2nd_sys[6];

  TFile *File_Input = TFile::Open("/Users/gavinwilks/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/output/RapidityDependentPaper/rho00y_all.root");
  for(int i_energy = 0; i_energy < 2; ++i_energy)
  {
    string GrapName_1st_stat = Form("g_rho00_order%d_%s_%s_StatError",1,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_1st_stat[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_1st_stat.c_str());
    string GrapName_1st_sys = Form("g_rho00_order%d_%s_%s_SysError",1,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_1st_sys[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_1st_sys.c_str());

    string GrapName_2nd_stat = Form("g_rho00_order%d_%s_%s_StatError",2,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_2nd_stat[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_2nd_stat.c_str());
    for(int i_point = 0; i_point < g_rho_2nd_stat[i_energy]->GetN(); ++i_point)
    {
      g_rho_2nd_stat[i_energy]->SetPointEXlow(i_point,0.0);
      g_rho_2nd_stat[i_energy]->SetPointEXhigh(i_point,0.0);
    }

    string GrapName_2nd_sys = Form("g_rho00_order%d_%s_%s_SysError",2,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_2nd_sys[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_2nd_sys.c_str());
  }


  TCanvas* c_rho00_1st = new TCanvas("c_rho00_1st","c_rho00_1st",10,20,450,960);
  c_rho00_1st->SetFillColor(10);
  c_rho00_1st->SetTopMargin(0.2);
  c_rho00_1st->SetBottomMargin(0.22);
  c_rho00_1st->SetRightMargin(0.05);
  c_rho00_1st->SetLeftMargin(0.15);
  c_rho00_1st->Divide(1,3,0.0,0.0,10);

  TH1F* h_frame_1st[3];
  for(int i_row = 0; i_row < 2; ++i_row)
  {
    string HistName = Form("h_frame_1st_%d",i_row);
    h_frame_1st[i_row] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.1,1.1);
    for(int i_bin = 0; i_bin < h_frame_1st[i_row]->GetNbinsX(); i_bin++)
    {
      h_frame_1st[i_row]->SetBinContent(i_bin+1,-10);
    }
  }

  for(int i_row = 0; i_row < 2; i_row++)
  {
    float scaling_factor = 1.0;
    //if(i_row == 2) scaling_factor = 0.78;
    c_rho00_1st->cd(i_row+1)->SetFillColor(10);
    c_rho00_1st->cd(i_row+1)->SetRightMargin(0.02);
    c_rho00_1st->cd(i_row+1)->SetLeftMargin(0.18);

    c_rho00_1st->cd(i_row+1);
    c_rho00_1st->cd(i_row+1)->SetTicks(1,1);
    c_rho00_1st->cd(i_row+1)->SetGrid(0,0);
    c_rho00_1st->cd(i_row+1)->SetFillColor(10);

    h_frame_1st[i_row]->SetStats(0);
    h_frame_1st[i_row]->SetTitle("");

    h_frame_1st[i_row]->GetXaxis()->SetTitle("|y|");
    h_frame_1st[i_row]->GetXaxis()->SetTitleSize(font_size*scaling_factor);
    h_frame_1st[i_row]->GetXaxis()->SetTitleOffset(0.85/scaling_factor);
    h_frame_1st[i_row]->GetXaxis()->CenterTitle();
    h_frame_1st[i_row]->GetXaxis()->SetRangeUser(pt_low,pt_high);
    h_frame_1st[i_row]->GetXaxis()->SetNdivisions(510,'N');
    h_frame_1st[i_row]->GetXaxis()->SetLabelSize(font_size*scaling_factor);

    h_frame_1st[i_row]->GetYaxis()->SetTitle("#rho_{00}");
    h_frame_1st[i_row]->GetYaxis()->CenterTitle();
    h_frame_1st[i_row]->GetYaxis()->SetTitleSize(font_size*scaling_factor);
    h_frame_1st[i_row]->GetYaxis()->SetTitleOffset(0.72/scaling_factor);
    h_frame_1st[i_row]->GetYaxis()->SetRangeUser(0.29,0.53);//rho00_low[i_row*2],rho00_high[i_row*2]);
    // h_frame_1st[i_row]->GetYaxis()->SetRangeUser(rho00_low[i_row],rho00_high[i_row]);
    h_frame_1st[i_row]->GetYaxis()->SetNdivisions(505,'N');
    h_frame_1st[i_row]->GetYaxis()->SetLabelSize(font_size*scaling_factor);
    h_frame_1st[i_row]->GetYaxis()->SetLabelOffset(0.01/scaling_factor);
    h_frame_1st[i_row]->DrawCopy("h");
    PlotLine(pt_low,pt_high,1.0/3.0,1.0/3.0,1,3,2);
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_1st_stat[i_row],plot_style_1st,plot_color_1st,0.0,plot_size_1st);
    plotSysErrors(g_rho_1st_sys[i_row],plot_color_1st);
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat[i_row],plot_style_2nd,plot_color_2nd,0.0,plot_size_2nd);
    plotSysErrors(g_rho_2nd_sys[i_row],plot_color_2nd);
    plotTopLegend((char*)mBeanEnergy[i_row].c_str(),3.6,0.28,font_size*scaling_factor,1,0.0,42,0);
    // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_1st_stat[i_row],plot_style_1st,plot_color_1st,plot_size_1st);
    // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat[i_row],plot_style_2nd,plot_color_2nd,plot_size_2nd);
    // plotTopLegend((char*)mBeanEnergy[i_row].c_str(),3.6,0.28,font_size*scaling_factor,1,0.0,42,0);

    if(i_row == 0)
    {
      plotTopLegend((char*)"Au+Au (0-80\% & |#eta| < 1)",1.5,0.48,leg_size*scaling_factor,1,0.0,42,0);

      Draw_TGAE_Point_new_Symbol(3.0,0.45,0.0,0.0,0.0,0.0,plot_style_1st,plot_color_1st,plot_size_1st);
      plotTopLegend((char*)"1^{st}-order EP",3.2,0.443,leg_size*scaling_factor,1,0.0,42,0);

      Draw_TGAE_Point_new_Symbol(3.0,0.42,0.0,0.0,0.0,0.0,plot_style_2nd,plot_color_2nd,plot_size_2nd);
      plotTopLegend((char*)"2^{nd}-order EP",3.2,0.413,leg_size*scaling_factor,1,0.0,42,0);
    }
  }

  c_rho00_1st->SaveAs("figures/c_rhoSys_y_BESII.pdf");

  //TCanvas* c_rho00_2nd = new TCanvas("c_rho00_2nd","c_rho00_2nd",500,20,450,960);
  //c_rho00_2nd->SetFillColor(10);
  //c_rho00_2nd->SetTopMargin(0.2);
  //c_rho00_2nd->SetBottomMargin(0.22);
  //c_rho00_2nd->SetRightMargin(0.05);
  //c_rho00_2nd->SetLeftMargin(0.15);
  //c_rho00_2nd->Divide(1,3,0.0,0.0,10);

  //TH1F* h_frame_2nd[3];
  //for(int i_row = 0; i_row < 3; ++i_row)
  //{
  //  string HistName = Form("h_frame_2nd_%d",i_row);
  //  h_frame_2nd[i_row] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,6);
  //  for(int i_bin = 0; i_bin < h_frame_2nd[i_row]->GetNbinsX(); i_bin++)
  //  {
  //    h_frame_2nd[i_row]->SetBinContent(i_bin+1,-10);
  //  }
  //}

  //for(int i_row = 0; i_row < 3; i_row++)
  //{
  //  float scaling_factor = 1.0;
  //  if(i_row == 2) scaling_factor = 0.78;
  //  c_rho00_2nd->cd(i_row+1)->SetFillColor(10);
  //  c_rho00_2nd->cd(i_row+1)->SetRightMargin(0.02);
  //  c_rho00_2nd->cd(i_row+1)->SetLeftMargin(0.18);

  //  c_rho00_2nd->cd(i_row+1);
  //  c_rho00_2nd->cd(i_row+1)->SetTicks(1,1);
  //  c_rho00_2nd->cd(i_row+1)->SetGrid(0,0);
  //  c_rho00_2nd->cd(i_row+1)->SetFillColor(10);

  //  h_frame_2nd[i_row]->SetStats(0);
  //  h_frame_2nd[i_row]->SetTitle("");

  //  h_frame_2nd[i_row]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //  h_frame_2nd[i_row]->GetXaxis()->SetTitleSize(font_size*scaling_factor);
  //  h_frame_2nd[i_row]->GetXaxis()->SetTitleOffset(0.85/scaling_factor);
  //  h_frame_2nd[i_row]->GetXaxis()->CenterTitle();
  //  h_frame_2nd[i_row]->GetXaxis()->SetRangeUser(pt_low,pt_high);
  //  h_frame_2nd[i_row]->GetXaxis()->SetNdivisions(510,'N');
  //  h_frame_2nd[i_row]->GetXaxis()->SetLabelSize(font_size*scaling_factor);

  //  h_frame_2nd[i_row]->GetYaxis()->SetTitle("#rho_{00}");
  //  h_frame_2nd[i_row]->GetYaxis()->CenterTitle();
  //  h_frame_2nd[i_row]->GetYaxis()->SetTitleSize(font_size*scaling_factor);
  //  h_frame_2nd[i_row]->GetYaxis()->SetTitleOffset(0.72/scaling_factor);
  //  h_frame_2nd[i_row]->GetYaxis()->SetRangeUser(rho00_low[i_row*2+1],rho00_high[i_row*2+1]);
  //  // h_frame_2nd[i_row]->GetYaxis()->SetRangeUser(rho00_low[i_row+3],rho00_high[i_row+3]);
  //  h_frame_2nd[i_row]->GetYaxis()->SetNdivisions(505,'N');
  //  h_frame_2nd[i_row]->GetYaxis()->SetLabelSize(font_size*scaling_factor);
  //  h_frame_2nd[i_row]->GetYaxis()->SetLabelOffset(0.01/scaling_factor);
  //  h_frame_2nd[i_row]->DrawCopy("h");
  //  PlotLine(pt_low,pt_high,1.0/3.0,1.0/3.0,1,3,2);
  //  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_1st_stat[i_row*2+1],plot_style_1st,plot_color_1st,plot_size_1st);
  //  plotSysErrors(g_rho_1st_sys[i_row*2+1],plot_color_1st);
  //  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat[i_row*2+1],plot_style_2nd,plot_color_2nd,plot_size_2nd);
  //  plotSysErrors(g_rho_2nd_sys[i_row*2+1],plot_color_2nd);
  //  plotTopLegend((char*)mBeanEnergy[i_row*2+1].c_str(),3.6,0.28,font_size*scaling_factor,1,0.0,42,0);
  //  // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_1st_stat[i_row+3],plot_style_1st,plot_color_1st,plot_size_1st);
  //  // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat[i_row+3],plot_style_2nd,plot_color_2nd,plot_size_2nd);
  //  // plotTopLegend((char*)mBeanEnergy[i_row+3].c_str(),3.6,0.28,font_size*scaling_factor,1,0.0,42,0);
  //}

  //c_rho00_2nd->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/c_rhoSys_pt_2nd.eps");
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    PlotLine(pt-0.1,pt+0.1,rho+err,rho+err,plot_color,2,1);
    PlotLine(pt-0.1,pt-0.1,rho+err-0.005,rho+err,plot_color,2,1);
    PlotLine(pt+0.1,pt+0.1,rho+err-0.005,rho+err,plot_color,2,1);
    PlotLine(pt-0.1,pt+0.1,rho-err,rho-err,plot_color,2,1);
    PlotLine(pt-0.1,pt-0.1,rho-err+0.005,rho-err,plot_color,2,1);
    PlotLine(pt+0.1,pt+0.1,rho-err+0.005,rho-err,plot_color,2,1);
  }
}
