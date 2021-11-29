#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/functions.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotSideBandCheck_27GeV()
{
  TCanvas *c_SideBand = new TCanvas("c_SideBand","c_SideBand",10,10,1000,500);
  c_SideBand->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_SideBand->cd(i_pad+1);
    c_SideBand->cd(i_pad+1)->SetLeftMargin(0.15);
    c_SideBand->cd(i_pad+1)->SetBottomMargin(0.15);
    c_SideBand->cd(i_pad+1)->SetTicks(1,1);
    c_SideBand->cd(i_pad+1)->SetGrid(0,0);
  }

  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/imass_27GeV.root");
  TH1D *h_mMass_SE = (TH1D*)File_Input->Get("imass_sig_before")->Clone("h_mMass_SE");
  TH1D *h_mMass_ME = (TH1D*)File_Input->Get("imass_bg")->Clone("h_mMass_ME");
  TH1D *h_mMass_SM = (TH1D*)File_Input->Get("imass_sig_after")->Clone("h_mMass_SM");

  c_SideBand->cd(1);
  h_mMass_SE->SetTitle("");
  h_mMass_SE->SetStats(0);
  h_mMass_SE->GetXaxis()->SetTitle("M(K^{+},K^{-})(GeV/c^{2})");
  h_mMass_SE->GetXaxis()->CenterTitle();
  h_mMass_SE->GetXaxis()->SetTitleSize(0.06);
  h_mMass_SE->GetXaxis()->SetTitleOffset(0.9);
  h_mMass_SE->GetXaxis()->SetLabelSize(0.04);
  h_mMass_SE->SetNdivisions(505,"X");

  h_mMass_SE->GetYaxis()->SetTitle("Yields");
  h_mMass_SE->GetYaxis()->CenterTitle();
  h_mMass_SE->GetYaxis()->SetTitleOffset(1.14);
  h_mMass_SE->GetYaxis()->SetTitleSize(0.06);
  h_mMass_SE->GetYaxis()->SetLabelSize(0.04);
  h_mMass_SE->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE->GetMaximum(),1.15*h_mMass_SE->GetMaximum());
  h_mMass_SE->SetNdivisions(505,"Y");
  h_mMass_SE->SetMarkerStyle(24);
  h_mMass_SE->SetMarkerColor(1);
  h_mMass_SE->SetMarkerSize(1.4);
  h_mMass_SE->DrawCopy("PE");
  // PlotLine(vmsa::InvMass_low[0],vmsa::InvMass_high[0],0.0,0.0,1,5,2);

  h_mMass_ME->SetMarkerStyle(1);
  // h_mMass_ME->SetLineColor(kAzure+2);
  // h_mMass_ME->SetFillColor(kAzure+2);
  h_mMass_ME->SetLineColor(kGray+2);
  h_mMass_ME->SetFillColor(kGray);
  h_mMass_ME->SetFillStyle(3003);
  h_mMass_ME->DrawCopy("h same");

  h_mMass_SM->SetMarkerStyle(20);
  h_mMass_SM->SetMarkerColor(2);
  h_mMass_SM->SetMarkerSize(1.2);
  h_mMass_SM->DrawCopy("pE same");

  TF1 *f_sig = new TF1("f_sig",PolyBreitWigner,vmsa::BW_Start[0],vmsa::BW_Stop[0],5);
  f_sig->SetParameter(0,vmsa::InvMass[0]);
  f_sig->SetParLimits(0,vmsa::InvMass[0]-1.5*vmsa::Width[0],vmsa::InvMass[0]+1.5*vmsa::Width[0]);
  f_sig->SetParameter(1,vmsa::Width[0]);
  f_sig->SetParLimits(1,0.004,0.070);
  f_sig->SetParameter(2,1.0);
  f_sig->SetParameter(3,-1.0);
  f_sig->SetParameter(4,1.0);
  f_sig->SetParameter(2,h_mMass_SM->GetMaximum()/f_sig->GetMaximum());
  f_sig->SetRange(vmsa::BW_Start[0],vmsa::BW_Stop[0]);
  h_mMass_SM->Fit(f_sig,"MNR");
  f_sig->SetLineColor(4);
  f_sig->SetLineStyle(1);
  f_sig->SetLineWidth(2);
  f_sig->DrawCopy("l same");

  TF1 *f_bg = new TF1("f_bg",Poly,vmsa::BW_Start[0],vmsa::BW_Stop[0],2);
  for(int i_par = 0; i_par < 2;++i_par)
  {
    f_bg->SetParameter(i_par,f_sig->GetParameter(i_par+3));
  }
  f_bg->SetLineColor(4);
  f_bg->SetLineStyle(2);
  f_bg->SetLineWidth(4);
  f_bg->DrawCopy("l same");

  // calculate Sig/Bkg 
  h_mMass_SM->Add(f_bg,-1.0);
  int inte_sig_start = h_mMass_SM->FindBin(1.01);
  int inte_sig_stop  = h_mMass_SM->FindBin(1.03);
  int inte_bkg_start = h_mMass_ME->FindBin(1.01);
  int inte_bkg_stop  = h_mMass_ME->FindBin(1.03);
  float Inte_SM = h_mMass_SM->Integral(inte_sig_start,inte_sig_stop);
  float Inte_ME = h_mMass_ME->Integral(inte_bkg_start,inte_bkg_stop);
  cout << "inte_sig_start = " << inte_sig_start << ", inte_sig_stop = " << inte_sig_stop << ", Inte_SE = " << Inte_SM << endl;
  cout << "inte_bkg_start = " << inte_bkg_start << ", inte_bkg_stop = " << inte_bkg_stop << ", Inte_ME = " << Inte_ME << endl;
  cout << "Sig/Bkg = " << Inte_SM/Inte_ME << endl;

  TLegend *leg_peak = new TLegend(0.5,0.4,0.8,0.55);
  leg_peak->SetBorderSize(0);
  leg_peak->SetFillColor(10);
  leg_peak->AddEntry(h_mMass_SE,"Same Event (Sig+Bg)","p");
  leg_peak->AddEntry(h_mMass_ME,"Mixed Event (Bg)","f");
  leg_peak->AddEntry(h_mMass_SM,"Signal","p");
  leg_peak->Draw("same");

  string leg_energy = Form("Au+Au %s & 20-60%%",vmsa::mBeamEnergy[3].c_str());
  plotTopLegend((char*)leg_energy.c_str(),0.18,0.84,0.04,1,0.0,42,1);
  plotTopLegend((char*)"1.2 < p_{T} < 1.8 GeV/c",0.20,0.79,0.04,1,0.0,42,1);
  string SBRatio = Form("Sig/Bkg = %1.2f",Inte_SM/Inte_ME);
  plotTopLegend((char*)SBRatio.c_str(),0.20,0.74,0.04,1,0.0,42,1);

  double invMass_27GeV = 1.02;
  double rhoData_27GeV = 0.360235;
  double errData_27GeV = 0.0025245586;
  TGraphErrors *g_rhoData_27GeV = new TGraphErrors();
  g_rhoData_27GeV->SetPoint(0,invMass_27GeV,rhoData_27GeV);
  g_rhoData_27GeV->SetPointError(0,0.01,errData_27GeV);

  double invMass[3] = {1.00,1.02,1.04};
  double errX[3] = {0.01,0.01,0.01};

  // double rhoData[3] = {0.347066,0.337032,0.317486};
  // double errData[3] = {0.00290068,0.00274057,0.00270384};
  double rhoData[3] = {0.344993,0.337498,0.319533};
  double errData[3] = {0.00173669,0.00168262,0.00166905};
  TGraphErrors *g_rhoData = new TGraphErrors();
  for(int i_point = 0; i_point < 3; ++i_point)
  {
    g_rhoData->SetPoint(i_point,invMass[i_point],rhoData[i_point]);
    g_rhoData->SetPointError(i_point,0.0,errData[i_point]);
  }

  double rhoBkg = 0.5*(rhoData[0]+rhoData[2]);
  double errBkg = 0.5*ErrorAdd(errData[0],errData[2]);
  TGraphErrors *g_rhoBkg = new TGraphErrors();
  g_rhoBkg->SetPoint(0,invMass[1],rhoBkg);
  g_rhoBkg->SetPointError(0,0.0,errBkg);

  double rhoSigBkg = (Inte_SM*rhoData_27GeV+Inte_ME*rhoBkg)/(Inte_SM+Inte_ME);
  double errSigBkg = ErrorAdd(Inte_SM*errData_27GeV,Inte_ME*errBkg)/(Inte_SM+Inte_ME); 
  TGraphErrors *g_rhoSigBkg = new TGraphErrors();
  g_rhoSigBkg->SetPoint(0,invMass[1]-0.002,rhoSigBkg);
  g_rhoSigBkg->SetPointError(0,0.0,errSigBkg);

  TH1F *h_frame = new TH1F("h_frame","h_frame",200,0.98,1.08);
  for(int i_bin = 0; i_bin < 200; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  c_SideBand->cd(2);
  h_frame->SetTitle("Au+Au 27GeV (Run11+Run18)");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.98,1.06);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("M(K^{+},K^{-})(GeV/c^{2})");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.29,0.37);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.98,1.06,1.0/3.0,1.0/3.0,1,2,2);

  g_rhoData_27GeV->SetMarkerStyle(29);
  g_rhoData_27GeV->SetMarkerColor(2);
  g_rhoData_27GeV->SetLineColor(2);
  g_rhoData_27GeV->SetMarkerSize(1.8);
  g_rhoData_27GeV->Draw("pE same");

  g_rhoData->SetMarkerStyle(20);
  g_rhoData->SetMarkerColor(kGray+2);
  g_rhoData->SetLineColor(kGray+2);
  g_rhoData->SetMarkerSize(1.4);
  g_rhoData->Draw("pE same");
  PlotLine(invMass[0],invMass[2],rhoData[0],rhoData[2],1,1,2);

  g_rhoBkg->SetMarkerStyle(30);
  g_rhoBkg->SetMarkerColor(4);
  g_rhoBkg->SetLineColor(4);
  g_rhoBkg->SetMarkerSize(1.4);
  g_rhoBkg->Draw("pE same");

  g_rhoSigBkg->SetMarkerStyle(30);
  g_rhoSigBkg->SetMarkerColor(2);
  g_rhoSigBkg->SetLineColor(2);
  g_rhoSigBkg->SetMarkerSize(1.4);
  g_rhoSigBkg->Draw("pE same");

  PlotLine(0.99,0.99,0.29,0.37,1,1,2);
  PlotLine(1.01,1.01,0.29,0.37,1,1,2);
  PlotLine(1.03,1.03,0.29,0.37,1,1,2);
  PlotLine(1.05,1.05,0.29,0.37,1,1,2);

  TLegend *leg = new TLegend(0.2,0.2,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(g_rhoData_27GeV,"STAR Signal","P");
  leg->AddEntry(g_rhoData,"#rho_{00} from different Inv. Mass","P");
  leg->AddEntry(g_rhoBkg,"Projected Bkg","P");
  leg->AddEntry(g_rhoSigBkg,"Sig & Bkg Averaged","P");
  leg->Draw("same");

  c_SideBand->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_SideBandCheck_27GeV.eps");
}
