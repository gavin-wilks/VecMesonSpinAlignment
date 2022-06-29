#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegend.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotSys_1stEnergyInte()
{
  const int mColor[3] = {1,4,6};
  const int mStyle[4] = {24,25,26,32};
  const string mMode[4] = {"Sigma_2_Inte","Sigma_0_Count","Sigma_1_Count","Sigma_2_Count"};
  const double pt_shift[3] = {-0.1,0.1,0.2};

  const int dca_default = 0;
  const double dca_center = 0.5;

  const int sig_default = 1;
  const double sig_center = 1.5;

  const int norm_default = 0;
  const double norm_center = 2.5;

  const int eff_default = 0;
  const double eff_center = 3.5;

  const double total_center = 4.5;
  // const int mode_default = 0;

  string inputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/NewF_JHChen/rho00_1stEnergyInte.root";
  cout << "Open InPut File: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  //==============================================================
  //--------------------------------------------------------------
  // get default value
  // const string GraphName_Default = "EP_1_eff_0_Dca_0_Sig_1_Phi_Norm_0_Sigma_2_Inte";
  const string GraphName_Default = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_2_Inte",eff_default,dca_default,sig_default,norm_default);
  TGraphAsymmErrors *g_rho_default = (TGraphAsymmErrors*)File_InPut->Get(GraphName_Default.c_str());
  cout << "Default Graph set to: " << GraphName_Default.c_str() << endl;

  double rho_def = 0.0;
  double err_stat = 0.0;
  TGraphAsymmErrors *g_rhoDef_Dca = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoDef_Sig = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoDef_Norm = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoDef_Eff = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoDef_Total = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoSys_Total = new TGraphAsymmErrors();
  for(int i_point = 0; i_point < 1; ++i_point)
  {
    double pt, rho;
    g_rho_default->GetPoint(i_point,pt,rho);
    g_rho_default->GetPoint(i_point,pt,rho_def);
    double err = g_rho_default->GetErrorYhigh(i_point);
    err_stat = g_rho_default->GetErrorYhigh(i_point);
    double width = 0.0;

    g_rhoDef_Dca->SetPoint(i_point,dca_center,rho);
    g_rhoDef_Dca->SetPointError(i_point,width,width,err,err);

    g_rhoDef_Sig->SetPoint(i_point,sig_center,rho);
    g_rhoDef_Sig->SetPointError(i_point,width,width,err,err);

    g_rhoDef_Norm->SetPoint(i_point,norm_center,rho);
    g_rhoDef_Norm->SetPointError(i_point,width,width,err,err);

    g_rhoDef_Eff->SetPoint(i_point,eff_center,rho);
    g_rhoDef_Eff->SetPointError(i_point,width,width,err,err);

    g_rhoDef_Total->SetPoint(i_point,total_center,rho);
    g_rhoDef_Total->SetPointError(i_point,width,width,err,err);

    g_rhoSys_Total->SetPoint(i_point,total_center,rho);
  }
  g_rhoDef_Dca->SetMarkerColor(2);
  g_rhoDef_Dca->SetMarkerStyle(29);
  g_rhoDef_Dca->SetMarkerSize(2.0);

  g_rhoDef_Sig->SetMarkerColor(2);
  g_rhoDef_Sig->SetMarkerStyle(29);
  g_rhoDef_Sig->SetMarkerSize(2.0);

  g_rhoDef_Norm->SetMarkerColor(2);
  g_rhoDef_Norm->SetMarkerStyle(29);
  g_rhoDef_Norm->SetMarkerSize(2.0);

  g_rhoDef_Eff->SetMarkerColor(2);
  g_rhoDef_Eff->SetMarkerStyle(29);
  g_rhoDef_Eff->SetMarkerSize(2.0);

  g_rhoDef_Total->SetMarkerColor(2);
  g_rhoDef_Total->SetMarkerStyle(29);
  g_rhoDef_Total->SetMarkerSize(2.0);
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // dca sysmatic value
  double rho_max_dca = rho_def;
  double rho_min_dca = rho_def;
  double rho_tmp_dca = 0.0;
  TGraphAsymmErrors *g_rhoSys_Dca[3][4];
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,i_dca,sig_default,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from DCA Cut: " << GraphName.c_str() << endl;

      TGraphAsymmErrors *g_rhoSys_Temp = (TGraphAsymmErrors*)File_InPut->Get(GraphName.c_str());
      g_rhoSys_Dca[i_dca][i_mode] = new TGraphAsymmErrors(); 
      for(int i_point = 0; i_point < 1; ++i_point) // 1st point is pT-integrated value
      {
	double pt, rho;
	g_rhoSys_Temp->GetPoint(i_point,pt,rho);
	double err = g_rhoSys_Temp->GetErrorYhigh(i_point);
	g_rhoSys_Dca[i_dca][i_mode]->SetPoint(i_point,dca_center+pt_shift[i_dca],rho);
	g_rhoSys_Dca[i_dca][i_mode]->SetPointError(i_point,0.0,0.0,err,err);
	rho_tmp_dca = rho;
      }
      g_rhoSys_Dca[i_dca][i_mode]->SetMarkerColor(mColor[i_dca]);
      g_rhoSys_Dca[i_dca][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      g_rhoSys_Dca[i_dca][i_mode]->SetMarkerSize(1.5);
      g_rhoSys_Dca[i_dca][i_mode]->SetLineColor(mColor[i_dca]);
      if(rho_tmp_dca > 0.0 && rho_tmp_dca < 1.0) 
      {
	if(rho_tmp_dca > rho_max_dca) rho_max_dca = rho_tmp_dca;
	if(rho_tmp_dca < rho_min_dca) rho_min_dca = rho_tmp_dca;
      }
    }
  }
  double delta_dca = rho_max_dca-rho_min_dca;
  double sys_dca = delta_dca/TMath::Sqrt(12.0);
  // string leg_dca = Form("#Delta#rho_{00,sys}^{dca}: %1.2f (#times 100%%)", 100.0*sys_dca/rho_def);
  string leg_dca = Form("#Delta#rho_{00,sys}^{dca}: %1.2f%%", 100.0*sys_dca/rho_def);
  cout << "sys_dca = " << sys_dca << endl;
  // cout << "sys_dca = " << sys_dca << " = " << 100.0*sys_dca/rho_def << "%" << endl;
  //--------------------------------------------------------------
  
  //--------------------------------------------------------------
  // sig sysmatic value
  double rho_max_sig = rho_def;
  double rho_min_sig = rho_def;
  double rho_tmp_sig = 0.0;
  TGraphAsymmErrors *g_rhoSys_Sig[3][4];
  for(int i_sig = 0; i_sig < 3; ++i_sig)
  {
    if(i_sig == sig_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,i_sig,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from nSigmaKaon Cut: " << GraphName.c_str() << endl;

      TGraphAsymmErrors *g_rhoSys_Temp = (TGraphAsymmErrors*)File_InPut->Get(GraphName.c_str());
      g_rhoSys_Sig[i_sig][i_mode] = new TGraphAsymmErrors();
      for(int i_point = 0; i_point < 1; ++i_point) // 1st point is pT-integrated value
      {
	double pt, rho;
	g_rhoSys_Temp->GetPoint(i_point,pt,rho);
	double err = g_rhoSys_Temp->GetErrorYhigh(i_point);
	g_rhoSys_Sig[i_sig][i_mode]->SetPoint(i_point,sig_center+pt_shift[i_sig],rho);
	g_rhoSys_Sig[i_sig][i_mode]->SetPointError(i_point,0.0,0.0,err,err);
	rho_tmp_sig = rho;
      }
      g_rhoSys_Sig[i_sig][i_mode]->SetMarkerColor(mColor[i_sig]);
      g_rhoSys_Sig[i_sig][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      g_rhoSys_Sig[i_sig][i_mode]->SetMarkerSize(1.5);
      g_rhoSys_Sig[i_sig][i_mode]->SetLineColor(mColor[i_sig]);
      if(rho_tmp_sig > 0.0 && rho_tmp_sig < 1.0) 
      {
	if(rho_tmp_sig > rho_max_sig) rho_max_sig = rho_tmp_sig;
	if(rho_tmp_sig < rho_min_sig) rho_min_sig = rho_tmp_sig;
      }
    }
  }
  double delta_sig = rho_max_sig - rho_min_sig;
  double sys_sig = delta_sig/TMath::Sqrt(12.0);
  cout << "sys_sig = " << sys_sig << endl;
  // string leg_sig = Form("#Delta#rho_{00,sys}^{sig}: %1.2f (#times 100%%)", 100.0*sys_sig/rho_def);
  string leg_sig = Form("#Delta#rho_{00,sys}^{sig}: %1.2f%%", 100.0*sys_sig/rho_def);
  // cout << "sys_sig = " << sys_sig << " = " << 100.0*sys_sig/rho_def << "%" << endl;
  //--------------------------------------------------------------
  
  //--------------------------------------------------------------
  // norm sysmatic value
  double rho_max_norm = rho_def;
  double rho_min_norm = rho_def;
  double rho_tmp_norm = 0.0;
  TGraphAsymmErrors *g_rhoSys_Norm[3][4];
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,sig_default,i_norm,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Normalization region: " << GraphName.c_str() << endl;

      TGraphAsymmErrors *g_rhoSys_Temp = (TGraphAsymmErrors*)File_InPut->Get(GraphName.c_str());
      g_rhoSys_Norm[i_norm][i_mode] = new TGraphAsymmErrors();
      for(int i_point = 0; i_point < 1; ++i_point) // 1st point is pT-integrated value
      {
	double pt, rho;
	g_rhoSys_Temp->GetPoint(i_point,pt,rho);
	double err = g_rhoSys_Temp->GetErrorYhigh(i_point);
	g_rhoSys_Norm[i_norm][i_mode]->SetPoint(i_point,norm_center+pt_shift[i_norm],rho);
	g_rhoSys_Norm[i_norm][i_mode]->SetPointError(i_point,0.0,0.0,err,err);
	rho_tmp_norm = rho;
      }
      g_rhoSys_Norm[i_norm][i_mode]->SetMarkerColor(mColor[i_norm]);
      g_rhoSys_Norm[i_norm][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      g_rhoSys_Norm[i_norm][i_mode]->SetMarkerSize(1.5);
      g_rhoSys_Norm[i_norm][i_mode]->SetLineColor(mColor[i_norm]);
      if(rho_tmp_norm > 0.0 && rho_tmp_norm < 1.0) 
      {
	if(rho_tmp_norm > rho_max_norm) rho_max_norm = rho_tmp_norm;
	if(rho_tmp_norm < rho_min_norm) rho_min_norm = rho_tmp_norm;
      }
    }
  }
  double delta_norm = rho_max_norm - rho_min_norm;
  double sys_norm = delta_norm/TMath::Sqrt(12.0);
  // string leg_norm = Form("#Delta#rho_{00,sys}^{norm}: %1.2f (#times 100%%)", 100.0*sys_norm/rho_def);
  string leg_norm = Form("#Delta#rho_{00,sys}^{norm}: %1.2f%%", 100.0*sys_norm/rho_def);
  cout << "sys_norm = " << sys_norm << endl;
  // cout << "sys_norm = " << sys_norm << " = " << 100.0*sys_norm/rho_def << "%" << endl;
  //--------------------------------------------------------------
  
  //--------------------------------------------------------------
  // eff sysmatic value
  double rho_max_eff = rho_def;
  double rho_min_eff = rho_def;
  double rho_tmp_eff = 0.0;
  TGraphAsymmErrors *g_rhoSys_Eff[2][4];
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    if(i_eff == eff_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",i_eff,dca_default,sig_default,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Eff: " << GraphName.c_str() << endl;

      TGraphAsymmErrors *g_rhoSys_Temp = (TGraphAsymmErrors*)File_InPut->Get(GraphName.c_str());
      g_rhoSys_Eff[i_eff][i_mode] = new TGraphAsymmErrors();
      for(int i_point = 0; i_point < 1; ++i_point) // 1st point is pT-integrated value
      {
	double pt, rho;
	g_rhoSys_Temp->GetPoint(i_point,pt,rho);
	double err = g_rhoSys_Temp->GetErrorYhigh(i_point);
	g_rhoSys_Eff[i_eff][i_mode]->SetPoint(i_point,eff_center+pt_shift[i_eff],rho);
	g_rhoSys_Eff[i_eff][i_mode]->SetPointError(i_point,0.0,0.0,err,err);
	rho_tmp_eff = rho;
      }
      g_rhoSys_Eff[i_eff][i_mode]->SetMarkerColor(mColor[i_eff]);
      g_rhoSys_Eff[i_eff][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      g_rhoSys_Eff[i_eff][i_mode]->SetMarkerSize(1.5);
      g_rhoSys_Eff[i_eff][i_mode]->SetLineColor(mColor[i_eff]);
      if(rho_tmp_eff > 0.0 && rho_tmp_eff < 1.0) 
      {
	if(rho_tmp_eff > rho_max_eff) rho_max_eff = rho_tmp_eff;
	if(rho_tmp_eff < rho_min_eff) rho_min_eff = rho_tmp_eff;
      }
    }
  }
  double delta_eff = rho_max_eff - rho_min_eff;
  double sys_eff = delta_eff/TMath::Sqrt(12.0);
  // string leg_eff = Form("#Delta#rho_{00,sys}^{eff}: %1.2f (#times 100%%)", 100.0*sys_eff/rho_def);
  string leg_eff = Form("#Delta#rho_{00,sys}^{eff}: %1.2f%%", 100.0*sys_eff/rho_def);
  cout << "sys_eff = " << sys_eff << endl;
  // cout << "sys_eff = " << sys_eff << " = " << 100.0*sys_eff/rho_def << "%" << endl;
  //--------------------------------------------------------------
  
  double sys_total = TMath::Sqrt(sys_dca*sys_dca + sys_sig*sys_sig + sys_norm*sys_norm + sys_eff*sys_eff);
  // string leg_total = Form("#Delta#rho_{00,sys}^{total}: %1.2f (#times 100%%)", 100.0*sys_total/rho_def);
  string leg_total = Form("#Delta#rho_{00,sys}^{total}: %1.2f%%", 100.0*sys_total/rho_def);
  cout << "sys_total = " << sys_total << endl;
  g_rhoSys_Total->SetPointError(0,0,0,sys_total,sys_total);
  // cout << "sys_total_graph = " << g_rhoSys_Total->GetErrorYhigh(0) << endl;

  cout << "rho_{00} = " << rho_def << " +/- " << err_stat << " +/- " << sys_total << endl;
  //==============================================================

  TCanvas *c_sys = new TCanvas("c_sys","c_sys",10,10,800,800);
  c_sys->cd();
  c_sys->cd()->SetLeftMargin(0.15);
  c_sys->cd()->SetBottomMargin(0.15);
  c_sys->cd()->SetTicks(1,1);
  c_sys->cd()->SetGrid(1,0);

  TH1F *h_frame = new TH1F("h_frame","h_frame",5,0.0,5.0);
  for(int bin_x = 1; bin_x < h_frame->GetNbinsX(); bin_x++)
  {
    h_frame->SetBinContent(bin_x,-10.0);
  }
  // string title_cuts = Form("Systematics @ AuAu %s",mBeamEnergy[energy].c_str());
  string title_cuts = "Systematics for Energy Integrated #rho_{00}";
  h_frame->SetTitle(title_cuts.c_str());
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetBinLabel(1,"dca");
  h_frame->GetXaxis()->SetBinLabel(2,"n#sigma_{K}");
  h_frame->GetXaxis()->SetBinLabel(3,"normalization");
  h_frame->GetXaxis()->SetBinLabel(4,"efficiency");
  h_frame->GetXaxis()->SetBinLabel(5,"total");
  // h_frame->GetXaxis()->LabelsOption("v");
  // h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  // h_frame->GetXaxis()->SetTitleSize(0.06);
  // h_frame->GetXaxis()->SetTitleOffset(1.1);
  // h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.33,0.40);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,3,2);

  // plot dca cuts
  g_rhoDef_Dca->Draw("pE same");
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Dca[i_dca][i_mode]->Draw("pE same");
    }
  }

  // plot sig cuts
  g_rhoDef_Sig->Draw("pE same");
  for(int i_sig = 0; i_sig < 3; ++i_sig)
  {
    if(i_sig == sig_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Sig[i_sig][i_mode]->Draw("pE same");
    }
  }

  // plot norm cuts
  g_rhoDef_Norm->Draw("pE same");
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Norm[i_norm][i_mode]->Draw("pE same");
    }
  }

  // plot eff cuts
  g_rhoDef_Eff->Draw("pE same");
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    if(i_eff == eff_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Eff[i_eff][i_mode]->Draw("pE same");
    }
  }

  // plot total sys errors
  g_rhoDef_Total->Draw("pE same");
  plotSysErrors(g_rhoSys_Total, 2);

  TLegend *leg = new TLegend(0.2,0.55,0.45,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(g_rhoDef_Total,"default","P");
  leg->AddEntry((TObject*)0,leg_dca.c_str(),"");
  leg->AddEntry((TObject*)0,leg_sig.c_str(),"");
  leg->AddEntry((TObject*)0,leg_norm.c_str(),"");
  leg->AddEntry((TObject*)0,leg_eff.c_str(),"");
  leg->AddEntry((TObject*)0,leg_total.c_str(),"");
  leg->Draw("same");

  string FigName = "/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/NewF_JHChen/c_SysTotal_1stEnergyInte.eps";
  c_sys->SaveAs(FigName.c_str());
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    PlotLine(pt-0.1,pt+0.1,rho+err,rho+err,plot_color,2,1);
    PlotLine(pt-0.1,pt-0.1,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(pt+0.1,pt+0.1,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(pt-0.1,pt+0.1,rho-err,rho-err,plot_color,2,1);
    PlotLine(pt-0.1,pt-0.1,rho-err+0.001,rho-err,plot_color,2,1);
    PlotLine(pt+0.1,pt+0.1,rho-err+0.001,rho-err,plot_color,2,1);
  }
}
