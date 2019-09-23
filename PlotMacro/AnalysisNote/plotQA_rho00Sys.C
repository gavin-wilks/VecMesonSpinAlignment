#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "../../Utility/draw.h"

void plotQA_rho00Sys(int energy = 6)
{
  string mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/rho00/RawRhoPtSys.root",mBeamEnergy[energy].c_str());
  cout << "inputfile is set to " << inputfile.c_str() << endl;

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_frame = (TH1F*)File_InPut->Get("h_frame");
  TGraphAsymmErrors *g_rho_count[3][3][3][3];
  TGraphAsymmErrors *g_rho_inte[3][3][3][3];

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  string leg_energy = Form("AuAu%s & 20-60 %%",mBeamEnergy[energy].c_str());
  h_frame->SetTitle(leg_energy.c_str());
  h_frame->GetYaxis()->SetRangeUser(0.25,0.4);
  h_frame->Draw("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    for(int i_sig = 0; i_sig < 3; ++i_sig)
    {
      for(int i_norm = 0; i_norm < 3; ++i_norm)
      {
	for(int i_sigma = 0; i_sigma < 3; ++i_sigma)
	{
	  string GraphName_Count = Form("rhoRaw_Centrality_9_2nd_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_%d_Count",i_dca,i_sig,i_norm,i_sigma);
	  // cout << "GraphName_Count = " << GraphName_Count.c_str() << endl;
	  g_rho_count[i_dca][i_sig][i_norm][i_sigma] = (TGraphAsymmErrors*)File_InPut->Get(GraphName_Count.c_str());
	  g_rho_count[i_dca][i_sig][i_norm][i_sigma]->SetMarkerStyle(24);
	  g_rho_count[i_dca][i_sig][i_norm][i_sigma]->SetMarkerColor(kGray+2);
	  g_rho_count[i_dca][i_sig][i_norm][i_sigma]->SetMarkerSize(1.2);
	  g_rho_count[i_dca][i_sig][i_norm][i_sigma]->Draw("pE same");

	  string GraphName_Inte = Form("rhoRaw_Centrality_9_2nd_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_%d_Inte",i_dca,i_sig,i_norm,i_sigma);
	  g_rho_inte[i_dca][i_sig][i_norm][i_sigma] = (TGraphAsymmErrors*)File_InPut->Get(GraphName_Inte.c_str());
	  g_rho_inte[i_dca][i_sig][i_norm][i_sigma]->SetMarkerStyle(24);
	  g_rho_inte[i_dca][i_sig][i_norm][i_sigma]->SetMarkerColor(2);
	  g_rho_inte[i_dca][i_sig][i_norm][i_sigma]->SetMarkerSize(1.2);
	  g_rho_inte[i_dca][i_sig][i_norm][i_sigma]->Draw("pE same");
	}
      }
    }
  }

  TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g_rho_count[0][0][0][0],"bin counting","p");
  leg->AddEntry(g_rho_inte[0][0][0][0],"bw integration","p");
  leg->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/RhoSig/c_rhoSys_AuAu%s.eps",mBeamEnergy[energy].c_str());
  c_rho->SaveAs(FigName.c_str());
}


