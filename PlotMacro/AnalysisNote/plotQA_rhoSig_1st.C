#include <TObject.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>
#include <TProfile2D.h>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"

void plotQA_rhoSig_1st(int energy = 6, int ptBin = 3)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AnalysisNote_1st/fig_%s.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();

  // read in histograms
  TH1F *h_mMass[7]; 
  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
  {
    string KEY = Form("imass_cos%d_pt%d",i_theta+1,ptBin);
    h_mMass[i_theta] = (TH1F*)File_InPut->Get(KEY.c_str());
  }
  string KEY_InteTheta = Form("imass_all_pt%d",ptBin);
  TH1F *h_mMass_InteTheta = (TH1F*)File_InPut->Get(KEY_InteTheta.c_str());

  string order[9] = {"I","II","III","IV","V","VI","VII","VIII","IX"};
  float mean = 0.0;
  float sigma = 0.0;
  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,900);
  c_diff->Divide(3,3);
  for(int i_theta = 0; i_theta < 9; ++i_theta)
  {
    c_diff->cd(i_theta+1);
    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
    c_diff->cd(i_theta+1)->SetTicks(1,1);
    c_diff->cd(i_theta+1)->SetGrid(0,0);
    if(i_theta < vmsa::CTS_stop)
    {
      string title = Form("%d/7 < cos(#theta*) < %d/7",i_theta,i_theta+1);
      h_mMass[i_theta]->SetTitle(title.c_str());
      h_mMass[i_theta]->SetStats(0);
      h_mMass[i_theta]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass[i_theta]->GetXaxis()->SetLabelSize(0.03);
      h_mMass[i_theta]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass[i_theta]->GetXaxis()->SetTitleSize(0.05);
      h_mMass[i_theta]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass[i_theta]->GetXaxis()->CenterTitle();

      h_mMass[i_theta]->GetYaxis()->SetRangeUser(h_mMass[i_theta]->GetMinimum(),1.1*h_mMass[i_theta]->GetMaximum());
      h_mMass[i_theta]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass[i_theta]->GetYaxis()->SetTitle("Yields");
      h_mMass[i_theta]->GetYaxis()->SetTitleSize(0.05);
      h_mMass[i_theta]->GetYaxis()->SetLabelSize(0.03);
      h_mMass[i_theta]->GetYaxis()->CenterTitle();

      h_mMass[i_theta]->SetMarkerStyle(24);
      h_mMass[i_theta]->SetMarkerColor(kGray+2);
      h_mMass[i_theta]->SetMarkerSize(1.2);
      h_mMass[i_theta]->Draw("pE");
      PlotLine(vmsa::InvMass_low[0],vmsa::InvMass_high[0],0.0,0.0,1,2,2);
      TLegend *leg_order;
      if(i_theta < 6)
      {
	leg_order = new TLegend(0.75,0.75,0.85,0.85);
      }
      else
      {
	leg_order = new TLegend(0.72,0.72,0.85,0.85);
      }
      leg_order->SetBorderSize(0);
      leg_order->SetFillColor(10);
      leg_order->AddEntry((TObject*)0,order[i_theta].c_str(),"");
      leg_order->Draw("same");
    }
    if(i_theta == vmsa::CTS_stop)
    {
      h_mMass_InteTheta->SetTitle("Integrated Yields");
      h_mMass_InteTheta->SetStats(0);
      h_mMass_InteTheta->GetXaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta->GetXaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass_InteTheta->GetXaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta->GetXaxis()->SetTitleOffset(1.2);
      h_mMass_InteTheta->GetXaxis()->CenterTitle();

      h_mMass_InteTheta->GetYaxis()->SetRangeUser(h_mMass_InteTheta->GetMinimum(),1.1*h_mMass_InteTheta->GetMaximum());
      h_mMass_InteTheta->GetYaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta->GetYaxis()->SetTitle("Yields");
      h_mMass_InteTheta->GetYaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta->GetYaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta->GetYaxis()->CenterTitle();

      h_mMass_InteTheta->SetMarkerStyle(24);
      h_mMass_InteTheta->SetMarkerColor(kGray+2);
      h_mMass_InteTheta->SetMarkerSize(1.2);
      h_mMass_InteTheta->Draw("pE");
      PlotLine(vmsa::InvMass_low[0],vmsa::InvMass_high[0],0.0,0.0,1,2,2);
      TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[0],vmsa::BW_Stop[0],3);
      f_bw->SetParameter(0,vmsa::InvMass[0]);
      f_bw->SetParLimits(0,vmsa::InvMass[0]-0.001,vmsa::InvMass[0]+0.001);
      f_bw->SetParameter(1,vmsa::Width[0]);
      f_bw->SetParameter(2,1.0);
      float norm = h_mMass_InteTheta->GetMaximum()/f_bw->GetMaximum();
      f_bw->SetParameter(2,norm);
      f_bw->SetRange(vmsa::BW_Start[0],vmsa::BW_Stop[0]);
      h_mMass_InteTheta->Fit(f_bw,"MQNR");
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(2);
      f_bw->SetLineWidth(2);
      f_bw->Draw("l same");
      mean = f_bw->GetParameter(0);
      sigma = f_bw->GetParameter(1);
      // cout << "mean = " << mean << ", sigma = " << sigma << endl;
      TLegend *leg_order = new TLegend(0.7,0.7,0.85,0.85);
      leg_order->SetBorderSize(0);
      leg_order->SetFillColor(10);
      leg_order->AddEntry((TObject*)0,order[i_theta].c_str(),"");
      leg_order->Draw("same");
    }
  }

  string KEY_Count = Form("fit_cos_pt%d",ptBin);
  TH1F *h_mYields_Count = (TH1F*)File_InPut->Get(KEY_Count.c_str());

  c_diff->cd(9);
  string leg_pt = Form("%1.1f < p_{T} < %1.1f GeV/c",vmsa::pt_low[energy][ptBin],vmsa::pt_up[energy][ptBin]);
  h_mYields_Count->SetTitle(leg_pt.c_str());
  h_mYields_Count->SetStats(0);
  h_mYields_Count->GetXaxis()->SetNdivisions(505,'N');
  h_mYields_Count->GetXaxis()->SetLabelSize(0.03);
  h_mYields_Count->GetXaxis()->SetTitle("cos(#theta*)");
  h_mYields_Count->GetXaxis()->SetTitleSize(0.05);
  h_mYields_Count->GetXaxis()->SetTitleOffset(1.2);
  h_mYields_Count->GetXaxis()->CenterTitle();

  h_mYields_Count->GetYaxis()->SetRangeUser(0.95*h_mYields_Count->GetMinimum(),1.05*h_mYields_Count->GetMaximum());
  h_mYields_Count->GetYaxis()->SetNdivisions(505,'N');
  h_mYields_Count->GetYaxis()->SetTitle("Counts");
  h_mYields_Count->GetYaxis()->SetTitleSize(0.05);
  h_mYields_Count->GetYaxis()->SetLabelSize(0.03);
  h_mYields_Count->GetYaxis()->CenterTitle();

  h_mYields_Count->SetMarkerStyle(24);
  h_mYields_Count->SetMarkerColor(4);
  h_mYields_Count->SetMarkerSize(1.2);
  h_mYields_Count->Draw("pE");

  TF1 *f_rho_count = new TF1("f_rho_count",SpinDensity,0.0,1.0,2);
  f_rho_count->SetParameter(0,0.33);
  f_rho_count->SetParameter(1,h_mYields_Count->GetMaximum());
  h_mYields_Count->Fit(f_rho_count,"NMRI");
  f_rho_count->SetLineColor(4);
  f_rho_count->SetLineWidth(2);
  f_rho_count->SetLineStyle(2);
  f_rho_count->Draw("l same");

  TLegend *leg = new TLegend(0.2,0.2,0.8,0.3);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(h_mYields_Count,"bin counting","p");
  leg->Draw("same");

  TLegend *leg_order = new TLegend(0.75,0.75,0.85,0.85);
  leg_order->SetBorderSize(0);
  leg_order->SetFillColor(10);
  leg_order->AddEntry((TObject*)0,order[8].c_str(),"");
  leg_order->Draw("same");

  string beamenergy = Form("AuAu %s", vmsa::mBeamEnergy[energy].c_str());
  plotTopLegend((char*)beamenergy.c_str(),0.2,0.78,0.07,1,0.0,42,1,1);
  plotTopLegend((char*)"1^{st} EP",0.3,0.68,0.07,1,0.0,42,1,1);

  for(int i_theta = 0; i_theta < 7; ++i_theta)
  {
    c_diff->cd(i_theta+1);
    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
    c_diff->cd(i_theta+1)->SetTicks(1,1);
    c_diff->cd(i_theta+1)->SetGrid(0,0);
    PlotLine(mean-2.0*sigma,mean-2.0*sigma,0.0,0.5*h_mMass[i_theta]->GetMaximum(),4,2,2);
    PlotLine(mean+2.0*sigma,mean+2.0*sigma,0.0,0.5*h_mMass[i_theta]->GetMaximum(),4,2,2);
  }

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/RhoSig/c_rhoSig_1st_AuAu%s_Pt%d.eps",vmsa::mBeamEnergy[energy].c_str(),ptBin);
  c_diff->SaveAs(FigName.c_str());
}
