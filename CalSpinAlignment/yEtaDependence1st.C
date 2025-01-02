#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void yEtaDependence1st(Int_t energy = 4, Int_t pid = 0, std::string setting = "AccRes")
{
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,10.05);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.2,1.2);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("|#eta| cut");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.32,0.4);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();


  TGraphAsymmErrors *g_before = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_after = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_beforeCent = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_afterCent = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_before1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_after1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_beforeCent1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_afterCent1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_vs = new TGraphAsymmErrors();


  float eta[4]       = {1.0,0.8,0.6,0.4};

  // Second Order Results
  float before[4]          = {0.3400,0.3441,0.3551,0.3817}; // 0-80%
  float beforeErr[4]       = {0.0008,0.0009,0.0010,0.0014}; // 0-80%
  float after[4]           = {0.3483,0.3438,0.3383,0.3366}; // v2 on and pt on, Eff + Acc, then Res correction 0-80%
  float afterErr[4]        = {0.0014,0.0015,0.0018,0.0024}; // v2 on and pt on, Eff + Acc, then Res correction 0-80%
  float beforeCent[4]    = {0.3399,0.3441,0.3545,0.3810}; // 0-80%
  float beforeErrCent[4] = {0.0008,0.0009,0.0011,0.0014}; // 0-80%
  float afterCent[4]     = {0.3494,0.3427,0.3370,0.3374}; // v2 on and pt on, Eff then Acc+Res correction 0-80%
  float afterErrCent[4]  = {0.0015,0.0016,0.0019,0.0025}; // v2 on and pt on, Eff then Acc+Res correction 0-80%

  // 1st Order Results
  float before1[4]        = {0.3376,0.3413,0.3518,0.3756}; // 0-80%
  float beforeErr1[4]     = {0.0007,0.0008,0.0009,0.0012}; // 0-80%
  float after1[4]         = {0.3467,0.3396,0.3326,0.3313}; // v2 on and pt on, Eff + Acc, then Res correction 0-80%
  float afterErr1[4]      = {0.0019,0.0021,0.0024,0.0031}; // v2 on and pt on, Eff + Acc, then Res correction 0-80%
  float beforeCent1[4]    = {0.3375,0.3411,0.3520,0.3778}; // 0-80%
  float beforeErrCent1[4] = {0.0007,0.0008,0.0009,0.0012}; // 0-80%
  float afterCent1[4]     = {0.3470,0.3386,0.3322,0.3324}; // v2 on and pt on, Eff then Acc+Res correction 0-80%
  float afterErrCent1[4]  = {0.0019,0.0026,0.0024,0.0031}; // v2 on and pt on, Eff then Acc+Res correction 0-80%


  for( int i = 0; i < 4; i++ )
  {
    g_after->SetPoint(i, eta[i], after[i]);
    g_after->SetPointError(i, 0.0, 0.0, afterErr[i], afterErr[i]);
    g_before->SetPoint(i, eta[i], before[i]);
    g_before->SetPointError(i, 0.0, 0.0, beforeErr[i], beforeErr[i]);
    g_afterCent->SetPoint(i, eta[i], afterCent[i]);
    g_afterCent->SetPointError(i, 0.0, 0.0, afterErrCent[i], afterErrCent[i]);
    g_beforeCent->SetPoint(i, eta[i], beforeCent[i]);
    g_beforeCent->SetPointError(i, 0.0, 0.0, beforeErrCent[i], beforeErrCent[i]);

    g_after1->SetPoint(i, eta[i], after1[i]);
    g_after1->SetPointError(i, 0.0, 0.0, afterErr1[i], afterErr1[i]);
    g_before1->SetPoint(i, eta[i], before1[i]);
    g_before1->SetPointError(i, 0.0, 0.0, beforeErr1[i], beforeErr1[i]);
    g_afterCent1->SetPoint(i, eta[i], afterCent1[i]);
    g_afterCent1->SetPointError(i, 0.0, 0.0, afterErrCent1[i], afterErrCent1[i]);
    g_beforeCent1->SetPoint(i, eta[i], beforeCent1[i]);
    g_beforeCent1->SetPointError(i, 0.0, 0.0, beforeErrCent1[i], beforeErrCent1[i]);

    g_vs->SetPoint(i, before[i], beforeCent[i]);
    g_vs->SetPointError(i, beforeErr[i], beforeErr[i], beforeErrCent[i], beforeErrCent[i]);
  }


  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  g_before->SetMarkerStyle(20);
  g_before->SetMarkerColor(kOrange+7);
  g_before->SetLineColor(kOrange+7);
  g_before->SetMarkerSize(1.3);
  g_before->Draw("pE same");
  
  g_after->SetMarkerStyle(20);
  g_after->SetMarkerColor(kBlue);
  g_after->SetLineColor(kBlue);
  g_after->SetMarkerSize(1.3);
  g_after->Draw("pE same");

  TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
  leg->AddEntry(g_before,"Before Corrections","p");
  leg->AddEntry(g_after,"After Corrections","p");
  leg->Draw("same");
  
  c_rho_SysError->SaveAs("figures/Phi/19GeV/yEtaDependence.pdf");   

  TCanvas *c_vs = new TCanvas("c_vs","c_vs",600,10,800,800);
  c_vs->cd();
  c_vs->cd()->SetLeftMargin(0.15);
  c_vs->cd()->SetBottomMargin(0.15);
  c_vs->cd()->SetTicks(1,1);
  c_vs->cd()->SetGrid(0,0);
  h_frame->SetTitle("After Corrections");
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  g_after->SetMarkerStyle(20);
  g_after->SetMarkerColor(kOrange+7);
  g_after->SetLineColor(kOrange+7);
  g_after->SetMarkerSize(1.3);
  g_after->Draw("pE same");

  g_afterCent->SetMarkerStyle(20);
  g_afterCent->SetMarkerColor(kBlue);
  g_afterCent->SetLineColor(kBlue);
  g_afterCent->SetMarkerSize(1.3);
  g_afterCent->Draw("pE same");

  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.8);
  leg2->AddEntry(g_after,"Rapidity Integrated Value","p");
  leg2->AddEntry(g_afterCent,"Centrality Integrated Value","p");
  leg2->Draw("same");

  c_vs->SaveAs("figures/Phi/19GeV/centEtaDependence_and_yEtaDependence_after.pdf");   


  TCanvas *c_vsb = new TCanvas("c_vsb","c_vsb",600,10,800,800);
  c_vsb->cd();
  c_vsb->cd()->SetLeftMargin(0.15);
  c_vsb->cd()->SetBottomMargin(0.15);
  c_vsb->cd()->SetTicks(1,1);
  c_vsb->cd()->SetGrid(0,0);
  h_frame->SetTitle("Before Corrections");
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  g_before->SetMarkerStyle(20);
  g_before->SetMarkerColor(kOrange+7);
  g_before->SetLineColor(kOrange+7);
  g_before->SetMarkerSize(1.3);
  g_before->Draw("pE same");

  g_beforeCent->SetMarkerStyle(20);
  g_beforeCent->SetMarkerColor(kBlue);
  g_beforeCent->SetLineColor(kBlue);
  g_beforeCent->SetMarkerSize(1.3);
  g_beforeCent->Draw("pE same");

  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.8);
  leg3->AddEntry(g_before,"Rapidity Integrated Value","p");
  leg3->AddEntry(g_beforeCent,"Centrality Integrated Value","p");
  leg3->Draw("same");

  c_vs->SaveAs("figures/Phi/19GeV/centEtaDependence_and_yEtaDependence_before.pdf");   

  TCanvas *c_vsao = new TCanvas("c_vsao","c_vsao",600,10,800,800);
  c_vsao->cd();
  c_vsao->cd()->SetLeftMargin(0.15);
  c_vsao->cd()->SetBottomMargin(0.15);
  c_vsao->cd()->SetTicks(1,1);
  c_vsao->cd()->SetGrid(0,0);
  h_frame->SetTitle("Before Corrections");
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  g_afterCent1->SetMarkerStyle(20);
  g_afterCent1->SetMarkerColor(kOrange+7);
  g_afterCent1->SetLineColor(kOrange+7);
  g_afterCent1->SetMarkerSize(1.3);
  g_afterCent1->Draw("pE same");

  g_afterCent->SetMarkerStyle(20);
  g_afterCent->SetMarkerColor(kBlue);
  g_afterCent->SetLineColor(kBlue);
  g_afterCent->SetMarkerSize(1.3);
  g_afterCent->Draw("pE same");

  TLegend *leg12 = new TLegend(0.5,0.6,0.8,0.8);
  leg12->AddEntry(g_afterCent1,"1^{st} Order Rapidity Integrated","p");
  leg12->AddEntry(g_afterCent, "2^{nd} Order Rapidity Integrated","p");
  leg12->Draw("same");

  c_vs->SaveAs("figures/Phi/19GeV/y_order_after.pdf");   
  
  
/*for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++) // Centrality loop
    {
      c_pt->cd(i_pt*3+i_cent+1);
      c_pt->cd(i_pt*3+i_cent+1)->SetLeftMargin(0.15);
      c_pt->cd(i_pt*3+i_cent+1)->SetBottomMargin(0.15);
      c_pt->cd(i_pt*3+i_cent+1)->SetTicks(1,1);
      c_pt->cd(i_pt*3+i_cent+1)->SetGrid(0,0);
      h_frame->DrawCopy("pE");
      PlotLine(-1.5,1.5,1.0/3.0,1.0/3.0,1,2,2);
      //Draw_TGAE_Point_new_Symbol(-0.1,0.46,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 

      //plotTopLegend((char*)leg_count.c_str(),0.0,0.457,0.03,1,0.0,42,0);

      //PlotLine(-0.25,-0.05,0.438,0.438,1,2,2);
      //plotTopLegend((char*)leg_line.c_str(),0.0,0.437,0.03,1,0.0,42,0);

      //plotTopLegend((char*)leg_energy.c_str(),-0.3,0.305,0.04,1,0.0,42,0);
      //plotTopLegend((char*)centStrings[i_cent].c_str(),-0.1,0.295,0.04,1,0.0,42,0);
      //plotTopLegend((char*)Form("%.2f<p_{T}<%.2f GeV/c",vmsa::pt_low_y[energy][i_pt],vmsa::pt_up_y[energy][i_pt]),-0.3,0.285,0.04,1,0.0,42,0);


      string KEY_Default = Form("rhoRaw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF0_Stat[KEY_Default],20,kBlue,0.0,1.0);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF1_Stat[KEY_Default],20,kOrange+7,0.0,1.0);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF2_Stat[KEY_Default],20,kGray+2,0.0,1.0);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF3_Stat[KEY_Default],20,kBlack,0.0,1.0);

      plotSysErrorsBox(g_mRhoF0_Sys[KEY_Default],kBlue);     
      plotSysErrorsBox(g_mRhoF1_Sys[KEY_Default],kOrange+7);     
      plotSysErrorsBox(g_mRhoF2_Sys[KEY_Default],kGray+2);     
      plotSysErrorsBox(g_mRhoF3_Sys[KEY_Default],kBlack);     

    }
  }
  c_pt->Update();
  c_pt->Print(outputname.c_str());
  string output_stop = Form("%s]",outputname.c_str());
  c_pt->Print(output_stop.c_str()); // close pdf file
  */

}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);
    
    //bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(vmsa::pt_low[beamE][i_energy],rho-err,vmsa::pt_up[beamE][i_energy],rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.08,rho-err,pt+0.08,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}
