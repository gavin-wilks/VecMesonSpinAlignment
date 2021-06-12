#include <string>

#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>

#include "../../../Utility/draw.h"
#include "../../../Utility/StSpinAlignmentCons.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotFig7_rho00PtPhi()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;

  const float size_marker = 1.4;
  
  float rho00_low[6] = {0.24,0.24,0.24,0.24,0.24,0.24}; // 11.5 -- 200 GeV
  float rho00_high[6] = {0.52,0.52,0.54,0.52,0.52,0.52};
  // float rho00_low[6] = {0.20,0.20,0.20,0.24,0.24,0.29}; // 11.5 -- 200 GeV
  // float rho00_high[6] = {0.50,0.52,0.64,0.40,0.40,0.35};
  float pt_low = 0.54;
  float pt_high = 5.54;
  float font_size = 0.07;
  float leg_size = 0.06;

  string mBeanEnergy[6] = {"11.5 GeV","19.6 GeV","27 GeV","39 GeV","62.4 GeV","200 GeV"};
  int mEnergy[6] = {11,19,27,39,62,200};

  TGraphAsymmErrors *g_rho_1st_stat_temp[6];
  TGraphAsymmErrors *g_rho_1st_sys_temp[6];
  TGraphAsymmErrors *g_rho_1st_stat[6];
  TGraphAsymmErrors *g_rho_1st_sys[6];
  TGraphAsymmErrors *g_rho_2nd_stat[6];
  TGraphAsymmErrors *g_rho_2nd_sys[6];

  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  for(int i_energy = 0; i_energy < 6; ++i_energy)
  {
    string GrapName_1st_stat = Form("rho00_1stEP_pt_stat_%d",mEnergy[i_energy]);
    g_rho_1st_stat_temp[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_1st_stat.c_str());
    g_rho_1st_stat[i_energy] = new TGraphAsymmErrors();
    for(int i_pt = 0; i_pt < g_rho_1st_stat_temp[i_energy]->GetN(); ++i_pt) // shift pT of 1st EP results by 0.1
    {
      double pt, rho;
      g_rho_1st_stat_temp[i_energy]->GetPoint(i_pt,pt,rho);
      double err = g_rho_1st_stat_temp[i_energy]->GetErrorYhigh(i_pt);
      g_rho_1st_stat[i_energy]->SetPoint(i_pt,pt+0.1,rho);
      g_rho_1st_stat[i_energy]->SetPointError(i_pt,0.0,0.0,err,err);
    }

    string GrapName_1st_sys = Form("rho00_1stEP_pt_sys_%d",mEnergy[i_energy]);
    g_rho_1st_sys_temp[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_1st_sys.c_str());
    g_rho_1st_sys[i_energy] = new TGraphAsymmErrors();
    for(int i_pt = 0; i_pt < g_rho_1st_sys_temp[i_energy]->GetN(); ++i_pt) // shift pT of 1st EP results by 0.1
    {
      double pt, rho;
      g_rho_1st_sys_temp[i_energy]->GetPoint(i_pt,pt,rho);
      double err = g_rho_1st_sys_temp[i_energy]->GetErrorYhigh(i_pt);
      g_rho_1st_sys[i_energy]->SetPoint(i_pt,pt+0.1,rho);
      g_rho_1st_sys[i_energy]->SetPointError(i_pt,0.0,0.0,err,err);
    }

    string GrapName_2nd_stat = Form("rho00_2ndEP_pt_stat_%d",mEnergy[i_energy]);
    g_rho_2nd_stat[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_2nd_stat.c_str());
    for(int i_point = 0; i_point < g_rho_2nd_stat[i_energy]->GetN(); ++i_point)
    {
      g_rho_2nd_stat[i_energy]->SetPointEXlow(i_point,0.0);
      g_rho_2nd_stat[i_energy]->SetPointEXhigh(i_point,0.0);
    }

    string GrapName_2nd_sys = Form("rho00_2ndEP_pt_sys_%d",mEnergy[i_energy]);
    g_rho_2nd_sys[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_2nd_sys.c_str());
  }

  //-----------------------------Set the Canvas---------------------------------
  const int N_x_pads = 2; //number of x-pads
  const int N_y_pads = 3; //number of y-pads
  const int N_total_pads = N_x_pads*N_y_pads; //number of total pads

  TCanvas *c_rho00_double = new TCanvas("c_rho00_double","c_rho00_double",10,10,800,1000);
  c_rho00_double->SetFillColor(10);
  c_rho00_double->Divide(N_x_pads,N_y_pads,0.0,0.0,10);

  //-----------------------------Set the coordinate-------------------------------  

  TH1F* h_frame[N_total_pads];
  for(int i_pad = 0; i_pad < N_total_pads; i_pad++)
  {
    string HistName = Form("h_frame_%d",i_pad);
    h_frame[i_pad] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,6.0);
    for(int bin_x = 1; bin_x < h_frame[i_pad]->GetNbinsX(); bin_x++)
    {
      h_frame[i_pad]->SetBinContent(bin_x,-10.0);
    }
    h_frame[i_pad]->GetXaxis()->SetRangeUser(0.5,5.5);
    h_frame[i_pad]->GetYaxis()->SetRangeUser(0.24,0.52);
    h_frame[i_pad]->SetStats(0);
    h_frame[i_pad]->SetTitle("");
    h_frame[i_pad]->GetXaxis()->SetNdivisions(510,'N');
    h_frame[i_pad]->GetYaxis()->SetNdivisions(505,'N');
    h_frame[i_pad]->GetYaxis()->SetLabelOffset(0.0);
    h_frame[i_pad]->GetXaxis()->SetLabelOffset(0.0);
    h_frame[i_pad]->GetYaxis()->SetLabelSize(0.0);
    h_frame[i_pad]->GetXaxis()->SetLabelSize(0.0);
  }

  //-----------------------------Set the coordinate--end--------------------------  

  double Margin_t = 0.01; // top offset in percent of total size
  double Margin_b = 0.1; // bottom offset in percent of total size
  double Margin_l = 0.19; // left offset in percent of total size
  double Margin_r = 0.01; // right offset in percent of total size
  double delta_y; // fraction in percent per pad in y direction
  double delta_x; // fraction in percent per pad in x direction
  double yUp_array[N_y_pads];
  double yDown_array[N_y_pads];
  double xLeft_array[N_x_pads];
  double xRight_array[N_x_pads];

  delta_y = (1.0-Margin_b-Margin_t)/(double)N_y_pads;
  delta_x = (1.0-Margin_l-Margin_r)/(double)N_x_pads;

  for(int y_pads = 0; y_pads < N_y_pads; y_pads++)//caculate pad size in y direction
  {
    yUp_array[y_pads] = 1.0 - Margin_t - y_pads*delta_y; 
    yDown_array[y_pads] = 1.0 - Margin_t - (y_pads+1.0)*delta_y; 
  }
  for(int x_pads = 0; x_pads < N_x_pads; x_pads++)
  {
    xLeft_array[x_pads] = Margin_l + x_pads*delta_x;
    xRight_array[x_pads] = Margin_l + (x_pads+1)*delta_x;
  }

  xLeft_array[0] = xLeft_array[0] - Margin_l;
  yDown_array[N_y_pads-1] = yDown_array[N_y_pads-1] - Margin_b;
  xRight_array[N_x_pads-1] = xRight_array[N_x_pads-1] + Margin_r;
  yUp_array[0] = yUp_array[0] + Margin_t;


  double scaling_pre = (yUp_array[N_y_pads-1]-yDown_array[N_y_pads-1])*(xRight_array[0]-xLeft_array[0]);
  //  double scaling_pre = (yUp_array[0]-yDown_array[0]);
  //  double scaling_pre = (xRight_array[0]-xLeft_array[0]);
  //    double scaling_pre = 1.0;

  for(int y_pads = 0; y_pads < N_y_pads; y_pads++)
  {
    for(int x_pads = 0; x_pads < N_x_pads; x_pads++)
    {
      int total_pad = (x_pads+1) + y_pads*N_x_pads;
      double scaling_factor = scaling_pre/((yUp_array[y_pads]-yDown_array[y_pads])*(xRight_array[x_pads]-xLeft_array[x_pads]));
      //      double scaling_factor = scaling_pre/(yUp_array[y_pads]-yDown_array[y_pads]);
      //      double scaling_factor = scaling_pre/(xRight_array[x_pads]-xLeft_array[x_pads]);
      //      cout << "scaling_factor = " << scaling_factor << endl;

      c_rho00_double->cd(total_pad);
      c_rho00_double->cd(total_pad)->SetTicks(1,1);
      c_rho00_double->cd(total_pad)->SetGrid(0,0);
      c_rho00_double->cd(total_pad)->SetFillColor(10);;
      c_rho00_double->cd(total_pad)->SetPad(xLeft_array[x_pads],yDown_array[y_pads],xRight_array[x_pads],yUp_array[y_pads]);


      if(x_pads == 0)
      {
	c_rho00_double->cd(total_pad)->SetLeftMargin(Margin_l/(xRight_array[x_pads]-xLeft_array[x_pads]));
	h_frame[total_pad-1]->GetYaxis()->SetLabelSize(0.08*scaling_factor);
	h_frame[total_pad-1]->GetYaxis()->SetLabelOffset(0.02*scaling_factor);
	h_frame[total_pad-1]->SetTickLength(0.02*scaling_factor,"X");
	h_frame[total_pad-1]->SetTickLength(0.02*scaling_factor,"Y");
      }
      if(x_pads == 0 && y_pads == N_y_pads-1)
      {
	c_rho00_double->cd(total_pad)->SetBottomMargin(Margin_b/(yUp_array[y_pads]-yDown_array[y_pads]));
	h_frame[total_pad-1]->SetTickLength(0.04*scaling_factor,"Y");
	h_frame[total_pad-1]->GetXaxis()->SetLabelSize(0.09*scaling_factor);
	h_frame[total_pad-1]->GetXaxis()->SetLabelOffset(0.01*scaling_factor);
	h_frame[total_pad-1]->GetYaxis()->SetLabelOffset(0.025*scaling_factor);
      }
      if(x_pads == N_x_pads-1)
      {
	c_rho00_double->cd(total_pad)->SetRightMargin(Margin_r/(xRight_array[x_pads]-xLeft_array[x_pads]));
	h_frame[total_pad-1]->GetXaxis()->SetLabelSize(0.065*scaling_factor);
	h_frame[total_pad-1]->GetXaxis()->SetLabelOffset(0.001*scaling_factor);
	h_frame[total_pad-1]->SetTickLength(0.01*scaling_factor,"X");
	h_frame[total_pad-1]->SetTickLength(0.02*scaling_factor,"Y");
      }
      if(x_pads == 1 && y_pads == N_y_pads-1)
      {
	c_rho00_double->cd(total_pad)->SetBottomMargin(Margin_b/(yUp_array[y_pads]-yDown_array[y_pads]));
	h_frame[total_pad-1]->SetTickLength(0.04*scaling_factor,"Y");
	h_frame[total_pad-1]->GetXaxis()->SetLabelSize(0.065*scaling_factor);
	h_frame[total_pad-1]->GetXaxis()->SetLabelOffset(0.003*scaling_factor);
      }
      if(y_pads == 0)
      {
	c_rho00_double->cd(total_pad)->SetTopMargin(Margin_t/(yUp_array[y_pads]-yDown_array[y_pads]));
      }

      if(x_pads == 0 && y_pads == 1)
      {
	// h_frame[total_pad-1]->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
	h_frame[total_pad-1]->GetYaxis()->SetTitle("#rho_{00}");
	h_frame[total_pad-1]->GetYaxis()->CenterTitle();
	h_frame[total_pad-1]->GetYaxis()->SetTitleSize(0.10*scaling_factor);
      }
      h_frame[total_pad-1]->DrawCopy("PE");
      PlotLine(pt_low,pt_high,1.0/3.0,1.0/3.0,1,3,2);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_1st_stat[total_pad-1],style_phi_1st,color_phi_1st,size_marker-0.2);
      plotSysErrors(g_rho_1st_sys[total_pad-1],color_phi_1st);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat[total_pad-1],style_phi_2nd,color_phi_2nd,size_marker+0.2);
      plotSysErrors(g_rho_2nd_sys[total_pad-1],color_phi_2nd);

      if(x_pads == 0 && y_pads != N_y_pads-1)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,font_size*scaling_factor,1,0.0,42,0);
      }
      if(x_pads == 0 && y_pads == N_y_pads-1)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,font_size*scaling_factor,1,0.0,42,0);
	h_frame[total_pad-1]->SetTickLength(0.03);
      }
      if(x_pads == N_x_pads-1 && y_pads == 1)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,0.050*scaling_factor,1,0.0,42,0);
      }
      if(x_pads == N_x_pads-1 && y_pads == N_y_pads-1)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,0.060*scaling_factor,1,0.0,42,0);
      }
      if(x_pads == N_x_pads-1 && y_pads == 0)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,0.05*scaling_factor,1,0.0,42,0);
      }
      if(x_pads == 0 && y_pads == 0)
      {
	plotTopLegend((char*)"Au+Au (20-60\% & |y| < 1.0)",0.9,0.49,leg_size*scaling_factor,1,0.0,42,0);

	Draw_TGAE_Point_new_Symbol(1.0,0.47,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
	plotTopLegend((char*)"#phi (1^{st}-order EP)",1.2,0.463,leg_size*scaling_factor,1,0.0,42,0);

	Draw_TGAE_Point_new_Symbol(1.0,0.44,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
	plotTopLegend((char*)"#phi (2^{nd}-order EP)",1.2,0.433,leg_size*scaling_factor,1,0.0,42,0);
      }

      // if(x_pads == 0 && y_pads == N_y_pads-1) plotTopLegend("#font[12]{p}_{#font[132]{T}}",0.94,0.10,0.1,1,0.0,42,1);
      // if(x_pads == 1 && y_pads == N_y_pads-1) plotTopLegend("#font[132]{(GeV}/#font[12]{c}#font[132]{)}",0.0,0.09,0.09,1,0.0,42,1);
      if(x_pads == 0 && y_pads == N_y_pads-1) plotTopLegend((char*)"p_{T}",0.94,0.10,0.1,1,0.0,42,1);
      if(x_pads == 1 && y_pads == N_y_pads-1) plotTopLegend((char*)"(GeV/c)",0.0,0.09,0.09,1,0.0,42,1);
    }
  }

  c_rho00_double->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig7_rho00PtPhi.eps");
  c_rho00_double->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig7_rho00PtPhi.png");
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
