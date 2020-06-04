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

void plotRho00InvMass()
{
  double sigma[6] = {2.0,2.5,3.0,1.5,1.0,0.5};
  double energy[6] = {11.0,19.6,27.0,39.0,62.4,200.0};
  double shift[6] = {0.0,-0.5,0.5,1.5,2.5,3.0};
  const int style[6] = {29,24,25,26,28,32};
  const int color[6] = {2,1,4,6,1,4};

  double rho_invmass0[6] = {0.385997,0.358477,0.362352,0.346422,0.355429,0.334372};
  double err_invmass0[6] = {0.0172128,0.00915756,0.00414334,0.00339564,0.00591969,0.00121709};

  double rho_invmass1[6] = {0.379131,0.359546,0.363309,0.347077,0.354366,0.334988};
  double err_invmass1[6] = {0.0181884,0.00965475,0.00413738,0.0035633,0.00623248,0.00121356};

  double rho_invmass2[6] = {0.38851,0.355959,0.368463,0.346431,0.355968,0.338578};
  double err_invmass2[6] = {0.0191525,0.010241,0.00414303,0.00371638,0.00653375,0.00121413};

  double rho_invmass3[6] = {0.382752,0.365186,0.366248,0.346136,0.356521,0.336299};
  double err_invmass3[6] = {0.0162102,0.00863566,0.00414508,0.00322725,0.00557174,0.00121813};

  double rho_invmass4[6] = {0.37722,0.364173,0.369286,0.345741,0.355271,0.336186};
  double err_invmass4[6] = {0.0153583,0.00815574,0.00414301,0.00307688,0.00531669,0.001215};

  double rho_invmass5[6] = {0.373272,0.368394,0.369168,0.343656,0.355511,0.336055};
  double err_invmass5[6] = {0.0156798,0.00836091,0.00413835,0.00312604,0.00544044,0.00121218};

  TGraphAsymmErrors *g_rhoInvMass[6];
  for(int i_sigma = 0; i_sigma < 6; ++i_sigma)
  {
    g_rhoInvMass[i_sigma] = new TGraphAsymmErrors();
    for(int i_energy = 0; i_energy < 6; ++i_energy)
    {
      if(i_sigma == 0) g_rhoInvMass[i_sigma]->SetPoint(i_energy,energy[i_energy]+shift[i_sigma],rho_invmass0[i_energy]);
      if(i_sigma == 0) g_rhoInvMass[i_sigma]->SetPointError(i_energy,0.0,0.0,err_invmass0[i_energy],err_invmass0[i_energy]);
      if(i_sigma == 1) g_rhoInvMass[i_sigma]->SetPoint(i_energy,energy[i_energy]+shift[i_sigma],rho_invmass1[i_energy]);
      if(i_sigma == 1) g_rhoInvMass[i_sigma]->SetPointError(i_energy,0.0,0.0,err_invmass1[i_energy],err_invmass0[i_energy]);
      if(i_sigma == 2) g_rhoInvMass[i_sigma]->SetPoint(i_energy,energy[i_energy]+shift[i_sigma],rho_invmass2[i_energy]);
      if(i_sigma == 2) g_rhoInvMass[i_sigma]->SetPointError(i_energy,0.0,0.0,err_invmass2[i_energy],err_invmass0[i_energy]);
      if(i_sigma == 3) g_rhoInvMass[i_sigma]->SetPoint(i_energy,energy[i_energy]+shift[i_sigma],rho_invmass3[i_energy]);
      if(i_sigma == 3) g_rhoInvMass[i_sigma]->SetPointError(i_energy,0.0,0.0,err_invmass3[i_energy],err_invmass0[i_energy]);
      if(i_sigma == 4) g_rhoInvMass[i_sigma]->SetPoint(i_energy,energy[i_energy]+shift[i_sigma],rho_invmass4[i_energy]);
      if(i_sigma == 4) g_rhoInvMass[i_sigma]->SetPointError(i_energy,0.0,0.0,err_invmass4[i_energy],err_invmass0[i_energy]);
      if(i_sigma == 5) g_rhoInvMass[i_sigma]->SetPoint(i_energy,energy[i_energy]+shift[i_sigma],rho_invmass5[i_energy]);
      if(i_sigma == 5) g_rhoInvMass[i_sigma]->SetPointError(i_energy,0.0,0.0,err_invmass5[i_energy],err_invmass0[i_energy]);
    }
  }

  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  c_rho00->cd()->SetLogx();
  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(9.0,240.0);
  h_frame->GetXaxis()->SetNdivisions(510,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.291,0.40);
  h_frame->GetYaxis()->SetNdivisions(510,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,1.0/3.0,1.0/3.0,1,3,2);

  for(int i_sigma = 0; i_sigma < 6; ++i_sigma)
  {
    g_rhoInvMass[i_sigma]->SetMarkerStyle(style[i_sigma]);
    g_rhoInvMass[i_sigma]->SetMarkerColor(color[i_sigma]);
    g_rhoInvMass[i_sigma]->SetMarkerSize(1.4);
    g_rhoInvMass[i_sigma]->SetLineColor(color[i_sigma]);
    g_rhoInvMass[i_sigma]->Draw("pE same");
  }

  double font_size = 0.035;
  plotTopLegend((char*)"Au+Au 20-60\%",13,0.392,font_size,1,0.0,42,0);
  plotTopLegend((char*)"1.2 < p_{T}< 5.4 GeV/c",37,0.392,font_size,1,0.0,42,0);
  plotTopLegend((char*)"|#eta| < 1",140,0.392,font_size,1,0.0,42,0);

  TLegend *leg = new TLegend(0.2,0.2,0.5,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  for(int i_sigma = 0; i_sigma < 6; ++i_sigma)
  {
    string invmass = Form("%1.1f#sigma", sigma[i_sigma]);
    leg->AddEntry(g_rhoInvMass[i_sigma],invmass.c_str(),"P");
  }
  leg->Draw("same");

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_rhoInvMass.eps");
}
