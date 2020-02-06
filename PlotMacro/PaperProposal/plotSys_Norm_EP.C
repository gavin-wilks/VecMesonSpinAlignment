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

using namespace std;

void plotSys_Norm_EP(int energy = 2)
{
  const string mBeamEnergy[6] = {"11.5 GeV","19.6 GeV","27 GeV","39 GeV","62.4 GeV","200 GeV"};
  const int mEnergy[6] = {11,19,27,39,62,200};
  const int mColor[3] = {1,4,6};
  const int mStyle[4] = {24,25,26,32};
  const string mMode[4] = {"Sigma_2_Inte","Sigma_0_Count","Sigma_1_Count","Sigma_2_Count"};
  const float pt_low = 0.54;
  const float pt_high = 5.54;
  const float pt_shift[3] = {-0.1,0.1,0.2};
  const int norm_default = 0;
  const string mLeg_norm[3] = {"[1.04,1.05]", "[0.99,1.0]", "both"};
  const string mLeg_mode[4] = {"BW Inte (2#sigma)", "Counting (2#sigma)", "Counting (2.5#sigma)", "Counting (3.0#sigma)"};

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/rho00_%dGeV_EP.root",mEnergy[energy]);
  if(energy == 2) inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/rho00_%dGeV_2ndMean_EP.root",mEnergy[energy]);
  cout << "Open InPut File: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  //--------------------------------------------------------------
  // get default value
  const string HistName_Default = "EP_2_eff_0_Dca_0_Sig_1_Phi_Norm_0_Sigma_2_Inte";
  TH1F *h_rho_default = (TH1F*)File_InPut->Get(HistName_Default.c_str());
  h_rho_default->SetMarkerColor(2);
  h_rho_default->SetMarkerStyle(29);
  h_rho_default->SetMarkerSize(2.0);
  cout << "Default Histogram set to: " << HistName_Default.c_str() << endl;

  TGraphAsymmErrors *g_rho_default = new TGraphAsymmErrors();
  for(int i_point = 0; i_point < h_rho_default->GetNbinsX(); ++i_point)
  {
    float pt = h_rho_default->GetBinCenter(i_point+1);
    float rho = h_rho_default->GetBinContent(i_point+1);
    float err = h_rho_default->GetBinError(i_point+1);
    float width = h_rho_default->GetBinWidth(i_point+1)/2.0;
    g_rho_default->SetPoint(i_point,pt,rho);
    g_rho_default->SetPointError(i_point,width,width,err,err);
  }
  g_rho_default->SetMarkerColor(2);
  g_rho_default->SetMarkerStyle(29);
  g_rho_default->SetMarkerSize(2.0);
  g_rho_default->RemovePoint(0); // 1st point is pT-integrated value
  //--------------------------------------------------------------

  TH1F *h_rhoSys_Norm[3][4]; // 0 for different Norm | 1 for different yields extraction
  TGraphAsymmErrors *g_rhoSys_Norm[3][4];
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_2_eff_0_Dca_0_Sig_1_Phi_Norm_%d_%s",i_norm,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from DCA Cut: " << HistName.c_str() << endl;
      h_rhoSys_Norm[i_norm][i_mode] = (TH1F*)File_InPut->Get(HistName.c_str());
      // h_rhoSys_Norm[i_norm][i_mode]->SetMarkerColor(mColor[i_norm]);
      h_rhoSys_Norm[i_norm][i_mode]->SetMarkerColor(1);
      h_rhoSys_Norm[i_norm][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      h_rhoSys_Norm[i_norm][i_mode]->SetMarkerSize(1.5);

      g_rhoSys_Norm[i_norm][i_mode] = new TGraphAsymmErrors();
      for(int i_point = 0; i_point < h_rho_default->GetNbinsX(); ++i_point)
      {
	float pt  = h_rhoSys_Norm[i_norm][i_mode]->GetBinCenter(i_point+1);
	float rho = h_rhoSys_Norm[i_norm][i_mode]->GetBinContent(i_point+1);
	float err = h_rhoSys_Norm[i_norm][i_mode]->GetBinError(i_point+1);
	g_rhoSys_Norm[i_norm][i_mode]->SetPoint(i_point,pt+pt_shift[i_norm],rho);
	g_rhoSys_Norm[i_norm][i_mode]->SetPointError(i_point,0.0,0.0,err,err);
      }
      g_rhoSys_Norm[i_norm][i_mode]->SetMarkerColor(mColor[i_norm]);
      g_rhoSys_Norm[i_norm][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      g_rhoSys_Norm[i_norm][i_mode]->SetMarkerSize(1.5);
      g_rhoSys_Norm[i_norm][i_mode]->SetLineColor(mColor[i_norm]);
      g_rhoSys_Norm[i_norm][i_mode]->RemovePoint(0);
    }
  }

  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,-0.5,6.0);
  for(int bin_x = 1; bin_x < h_frame->GetNbinsX(); bin_x++)
  {
    h_frame->SetBinContent(bin_x,-10.0);
  }
  string title_cuts = Form("Normalization Systematics @ AuAu %s",mBeamEnergy[energy].c_str());
  h_frame->SetTitle(title_cuts.c_str());
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(pt_low,pt_high);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.22,0.38);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (In-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(pt_low,pt_high,1.0/3.0,1.0/3.0,1,3,2);

  // h_rho_default->Draw("pE same");
  g_rho_default->Draw("pE same");
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      // h_rhoSys_Norm[i_norm][i_mode]->Draw("pE same");
      g_rhoSys_Norm[i_norm][i_mode]->Draw("pE same");
    }
  }

  string formula = "#Delta#rho_{00,sys}^{norm} = #frac{#rho_{00,max}^{norm}-#rho_{00,min}^{norm}}{#sqrt{12}}";
  // string formula = "#Delta#rho_{00,sys}^{norm} = (#rho_{00,max}^{norm}-#rho_{00,min}^{norm})/#sqrt{12}";
  plotTopLegend((char*)formula.c_str(),0.3,0.45,0.03,1,0.0,42,1);

  TLegend *leg_norm = new TLegend(0.2,0.2,0.5,0.4);
  leg_norm->SetBorderSize(0);
  leg_norm->SetFillColor(10);
  leg_norm->AddEntry(g_rho_default,"default ([1.04,1.05])","P");
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    leg_norm->AddEntry(g_rhoSys_Norm[i_norm][0],mLeg_norm[i_norm].c_str(),"P");
  }
  leg_norm->AddEntry((TObject*)0," ","");
  leg_norm->Draw("same");

  TLegend *leg_mode = new TLegend(0.5,0.2,0.7,0.4);
  leg_mode->SetBorderSize(0);
  leg_mode->SetFillColor(10);
  for(int i_mode = 0; i_mode< 4; ++i_mode)
  {
    leg_mode->AddEntry(h_rhoSys_Norm[1][i_mode],mLeg_mode[i_mode].c_str(),"P");
  }
  leg_mode->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_SysNorm_AuAu%dGeV_EP.eps",mEnergy[energy]);
  c_rho00->SaveAs(FigName.c_str());
}
