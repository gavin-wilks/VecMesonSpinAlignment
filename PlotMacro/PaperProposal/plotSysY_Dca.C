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

void plotSysY_Dca(int energy = 0, int order = 1, string correction = "AccRes", string etamode = "eta1_eta1")
{
  const string EP[2] = {"","2nd"};
  const string mBeamEnergy[2] = {"14.6 GeV","19.6 GeV"};
  const string mBeamEnergyFile[2] = {"14GeV","19GeV"};
  const int mEnergy[2] = {14,19};
  const int mColor[3] = {1,4,6};
  const int mStyle[4] = {24,25,26,32};
  const string mMode[4] = {"Sigma_0_Inte","Sigma_0_Count","Sigma_1_Count","Sigma_2_Count"};
  const float pt_low = 0.0;
  const float pt_high = 1.1;
  const float pt_shift[3] = {0,0.03,0.06};
  const int dca_default = 0;
  const string mLeg_dca[3] = {"dca < 2.0 cm", "dca < 2.5 cm", "dca < 3.0 cm"};
  const string mLeg_mode[4] = {"BW Inte (2#sigma)", "Counting (2#sigma)", "Counting (2.5#sigma)", "Counting (3.0#sigma)"};

  string inputfile = Form("../../output/AuAu%s/Phi/RhoEta_%sSysErrors_%s_PolySys.root",mBeamEnergyFile[energy].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../../output/AuAu%s/Phi/RhoEta_%sSysErrors_%s_PolySys_FirstOrder.root",mBeamEnergyFile[energy].c_str(),correction.c_str(),etamode.c_str());

  cout << "Open InPut File: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->Print();

  //--------------------------------------------------------------
  // get default value
  string HistName_Default = Form("rhoRaw_%s_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1",EP[order-1].c_str());
  TGraphAsymmErrors *g_rho_default = (TGraphAsymmErrors*)File_InPut->Get(HistName_Default.c_str());
  cout << "Loaded Default" << endl;
  cout << HistName_Default << endl;
  g_rho_default->Print();
  g_rho_default->SetMarkerColor(2);
  g_rho_default->SetMarkerStyle(29);
  g_rho_default->SetMarkerSize(2.0);
  cout << "SetMarker" << endl;
  //g_rho_default->RemovePoint(0); // 1st point is pT-integrated value
  //--------------------------------------------------------------
  
  cout << "Before Loading Different yield extractions " << endl;

  TGraphAsymmErrors *g_rhoSys_Dca[3][4]; // 0 for different Dca | 1 for different yields extraction
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      cout << "Before creating histname" << endl;
      string HistName = Form("rhoRaw_%s_Dca_%d_Sig_0_Phi_Norm_0_%s_Poly1",EP[order-1].c_str(),i_dca,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from DCA Cut: " << HistName.c_str() << endl;
      g_rhoSys_Dca[i_dca][i_mode] = (TGraphAsymmErrors*)File_InPut->Get(HistName.c_str());
      g_rhoSys_Dca[i_dca][i_mode]->SetMarkerColor(mColor[i_dca]);
      g_rhoSys_Dca[i_dca][i_mode]->SetMarkerStyle(mStyle[i_mode]);
      g_rhoSys_Dca[i_dca][i_mode]->SetMarkerSize(1.5);
      g_rhoSys_Dca[i_dca][i_mode]->SetLineColor(mColor[i_dca]);
      for(int i = 0; i < g_rhoSys_Dca[i_dca][i_mode]->GetN(); i++)
      {
        double x, y;
        g_rhoSys_Dca[i_dca][i_mode]->GetPoint(i,x,y);
        g_rhoSys_Dca[i_dca][i_mode]->SetPoint(i,x+pt_shift[i_dca],y);
      }
      //g_rhoSys_Dca[i_dca][i_mode]->RemovePoint(0);
    }
  }

  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,-0.5,100);
  for(int bin_x = 1; bin_x < h_frame->GetNbinsX(); bin_x++)
  {
    h_frame->SetBinContent(bin_x,-10.0);
  }
  string title_cuts = Form("DCA Cuts Systematics @ AuAu %s EP Order %d",mBeamEnergy[energy].c_str(),order);
  h_frame->SetTitle(title_cuts.c_str());
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(pt_low,pt_high);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("|y|");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.2,0.55);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(pt_low,pt_high,1.0/3.0,1.0/3.0,1,3,2);

  // h_rho_default->Draw("pE same");
  g_rho_default->Draw("pE same");
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      // h_rhoSys_Dca[i_dca][i_mode]->Draw("pE same");
      g_rhoSys_Dca[i_dca][i_mode]->Draw("pE same");
    }
  }

  string formula = "#Delta#rho_{00,sys}^{dca} = #frac{#rho_{00,max}^{dca}-#rho_{00,min}^{dca}}{#sqrt{12}}";
  // string formula = "#Delta#rho_{00,sys}^{dca} = (#rho_{00,max}^{dca}-#rho_{00,min}^{dca})/#sqrt{12}";
  plotTopLegend((char*)formula.c_str(),0.3,0.8,0.03,1,0.0,42,1);

  TLegend *leg_dca = new TLegend(0.2,0.2,0.5,0.4);
  leg_dca->SetBorderSize(0);
  leg_dca->SetFillColor(10);
  leg_dca->AddEntry(g_rho_default,"default (dca < 2.0 cm)","P");
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    leg_dca->AddEntry(g_rhoSys_Dca[i_dca][0],mLeg_dca[i_dca].c_str(),"P");
  }
  leg_dca->AddEntry((TObject*)0," ","");
  leg_dca->Draw("same");

  TLegend *leg_mode = new TLegend(0.5,0.2,0.7,0.4);
  leg_mode->SetBorderSize(0);
  leg_mode->SetFillColor(10);
  for(int i_mode = 0; i_mode< 4; ++i_mode)
  {
    leg_mode->AddEntry(g_rhoSys_Dca[1][i_mode],mLeg_mode[i_mode].c_str(),"P");
  }
  leg_mode->Draw("same");

  string FigName = Form("c_SysYDca_AuAu%dGeV_order%d.pdf",mEnergy[energy],order);
  c_rho00->SaveAs(FigName.c_str());
}
