#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compRhoEff_PhiOrderY_Folded(Int_t energy = 0, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw")
{
  double rapidity[11] = {0.,	                
                         0.1,	                
                         0.2,	                
                         0.30000000000000004,
                         0.4,	                
                         0.5,	                
                         0.6000000000000001,	
                         0.7000000000000001,	
                         0.8,	                
                         0.9,	                
                         1.};	                
  
  double rhoTheory[11] = {0.2979499013864536,
                          0.299443815399517,
                          0.3039886168598637,
                          0.3117758225131158,
                          0.3231325414288237,
                          0.3385337619154553,
                          0.35862038503066823,
                          0.38422384087303474,
                          0.41639836866905044,
                          0.45646236573878285,
                          0.5060504797489073};

  double rhoTheoryErr[11] = {0.026776,
	                     0.0266206,
	                     0.0262643,
	                     0.0260758,
	                     0.0267701,
	                     0.0293957,
                             0.035037,
	                     0.0444524,
	                     0.0581137,
	                     0.0764945,
                             0.100262};

  TGraphErrors *g_rhoTheory = new TGraphErrors();
  for(int i = 0; i < 11; i++)
  {
    g_rhoTheory->SetPoint(i,rapidity[i],rhoTheory[i]);
    g_rhoTheory->SetPointError(i,0.0,rhoTheoryErr[i]);
  }

  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  cout << "Before loading plots" << endl;

  string inputfile1 = Form("../output/AuAu%s/%s/RhoEta_%sSysErrors_eta1_eta1_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  string inputfile2 = Form("../output/AuAu%s/%s/RhoEta_%sSysErrors_eta1_eta1_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());

  TFile *File_InPut1 = TFile::Open(inputfile1.c_str());
  TFile *File_InPut2 = TFile::Open(inputfile2.c_str());

  cout << "Before loading plots for stat and sys" << endl;
  TString StatErrorRho1 = Form("g_rho00_order1_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString SysErrorRho1 = Form("g_rho00_order1_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString StatErrorRho2 = Form("g_rho00_order2_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString SysErrorRho2 = Form("g_rho00_order2_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TGraphAsymmErrors *g_Stat1   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPut1->Get(StatErrorRho1.Data()))->Clone();
  TGraphAsymmErrors *g_Sys1    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPut1->Get( SysErrorRho1.Data()))->Clone();
  TGraphAsymmErrors *g_Stat2   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPut2->Get(StatErrorRho2.Data()))->Clone();
  TGraphAsymmErrors *g_Sys2    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPut2->Get( SysErrorRho2.Data()))->Clone();

  cout << "Before loading h_frame" << endl;
  TH1F *h_frame = (TH1F*)File_InPut1->Get("h_frame");

  cout << "Loaded plots" << endl;

  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  //TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //cout << "All good" << endl;

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->SetTitle("");
  h_frame->GetXaxis()->SetRangeUser(0.0,1.0);
  h_frame->GetYaxis()->SetRangeUser(0.29,0.53);
  h_frame->GetXaxis()->SetTitle("|y|");
  //if(energy == 3) h_frame->GetYaxis()->SetRangeUser(0.20,0.56);
  h_frame->Draw("pE");
  //g_StatErrors->SetMarkerStyle(20);
  //g_StatErrors->SetMarkerColor(kGray+2);
  //g_StatErrors->SetLineColor(2);
  //g_StatErrors->Draw("pE same");
  //g_SysErrors->SetMarkerStyle(20);
  //g_SysErrors->SetMarkerColor(kGray+2);
  //g_SysErrors->SetLineColor(kGray+2);
  //g_SysErrors->Draw("pE same");
//  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  PlotLine(0.0,1.0,1.0/3.0,1.0/3.0,1,2,2);

  //besi19->SetLineColor(kBlack);
  //besi19->SetMarkerStyle(20);
  //besi19->SetMarkerSize(1.3);
  //besi19->SetMarkerColor(kBlack);
  //besi19->SetLineColor(kBlack);
  //besi19->Draw("pE same");

  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->SetMarkerStyle(20);
  //besi19_sys->SetMarkerSize(1.3);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //plotSysErrorsBox(besi19_sys,kBlack);//,energy);
  //besi19_sys->Draw("pE [] same");

  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  //TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //besi19->SetMarkerStyle(20);
  //besi19->SetMarkerSize(1.3);
  //besi19->SetMarkerColor(kBlack);
  //besi19->SetLineColor(kBlack);
  //besi19->Draw("pE Z same");

  //besi19_sys->SetMarkerStyle(20);
  //besi19_sys->SetMarkerSize(1.3);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->Draw("pE [] same");

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.49,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.487,0.03,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(0.5,0.47,0.0,0.0,0.0,0.0,20,kBlue,1.3); 
  //string leg_count0 = Form("#phi F");
  //plotTopLegend((char*)leg_count0.c_str(),0.6,0.467,0.03,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  //string leg_count1 = Form("#phi F+#deltaF");
  //plotTopLegend((char*)leg_count1.c_str(),0.6,0.447,0.03,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,kGray+2,1.3); 
  //string leg_count2 = Form("#phi F-#deltaF");
  //plotTopLegend((char*)leg_count2.c_str(),0.6,0.427,0.03,1,0.0,42,0);
  //Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  //string leg_countA = "Eff+Acc+Res Corrected #phi BES-II";
  //plotTopLegend((char*)leg_countA.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);
  if(energy == 4)
  {
    g_rhoTheory->SetLineWidth(5.0);
    g_rhoTheory->SetLineColor(kBlue);
    g_rhoTheory->SetLineStyle(2);
    g_rhoTheory->SetFillColorAlpha(kBlue,0.25);
    g_rhoTheory->SetFillStyle(3001);
    //g_rhoTheory->Draw("same C3");
    g_rhoTheory->Draw("same CX");
  }
  //if(energy == 4)
  {
    Draw_TGAE_Point_new_Symbol(0.1,0.50,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 
    string leg_count = "#phi 2^{nd} order (1.0 < p_{T} < 5.0 GeV/c)";
    plotTopLegend((char*)leg_count.c_str(),0.15,0.497,0.03,1,0.0,42,0);

    string leg_besi = "#phi 1^{st} order (1.0 < p_{T} < 5.0 GeV/c)";
    Draw_TGAE_Point_new_Symbol(0.1,0.485,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
    plotTopLegend((char*)leg_besi.c_str(),0.15,0.482,0.03,1,0.0,42,0);
 
    //string leg_sp = "STAR Preliminary";
    //plotTopLegend((char*)leg_sp.c_str(),0.1,0.515,0.03,2,0.0,42,0);

    PlotLine(0.05,0.125,0.47,0.47,1,2,2);
    string leg_line = "#rho_{00} = 1/3";
    plotTopLegend((char*)leg_line.c_str(),0.15,0.47,0.03,1,0.0,42,0);
    if(energy == 4)
    {
      PlotLine(0.05,0.125,0.455,0.455,kBlue,3,2);
      string leg_line2 = "Theory (1.2 < p_{T} < 5.4 GeV/c,";
      plotTopLegend((char*)leg_line2.c_str(),0.15,0.455,0.03,1,0.0,42,0);
      string leg_line2b = "              20%-60%)";
      plotTopLegend((char*)leg_line2b.c_str(),0.15,0.4425,0.03,1,0.0,42,0);
    }
    string leg_energy = Form("Au+Au %s", vmsa::mBeamEnergyText[energy].c_str());
    plotTopLegend((char*)leg_energy.c_str(),0.1,0.40,0.04,1,0.0,42,0);
    plotTopLegend((char*)"0%-80%",0.2,0.385,0.04,1,0.0,42,0);
  }
  //if(energy == 3)
  //{
  //  Draw_TGAE_Point_new_Symbol(-0.6,0.50,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 
  //  string leg_count = "#phi 2^{nd} order (1.0 < p_{T} < 5.0 GeV/c)";
  //  plotTopLegend((char*)leg_count.c_str(),-0.5,0.497,0.03,1,0.0,42,0);

  //  string leg_besi = "#phi 1^{st} order (1.0 < p_{T} < 5.0 GeV/c)";
  //  Draw_TGAE_Point_new_Symbol(-0.6,0.485,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  //  plotTopLegend((char*)leg_besi.c_str(),-0.5,0.482,0.03,1,0.0,42,0);
 
  //  string leg_sp = "STAR Preliminary";
  //  plotTopLegend((char*)leg_sp.c_str(),-0.6,0.515,0.03,2,0.0,42,0);

  //  PlotLine(-0.7,-0.55,0.47,0.47,1,2,2);
  //  string leg_line = "#rho_{00} = 1/3";
  //  plotTopLegend((char*)leg_line.c_str(),-0.5,0.47,0.03,1,0.0,42,0);

  //  string leg_energy = Form("Au+Au %s", vmsa::mBeamEnergyText[energy].c_str());
  //  plotTopLegend((char*)leg_energy.c_str(),-0.35,0.43,0.04,1,0.0,42,0);
  //  plotTopLegend((char*)"0%-80%",-0.15,0.415,0.04,1,0.0,42,0);
  //}

  //PlotLine(0.25,0.75,0.262,0.262,1,2,2);
  //string leg_line = "#rho_{00} = 1/3";
  //plotTopLegend((char*)leg_line.c_str(),0.85,0.26,0.03,1,0.0,42,0);


  //string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  //plotTopLegend((char*)leg_energy.c_str(),2.5,0.26,0.04,1,0.0,42,0);
  //plotTopLegend((char*)"20%-60%",3.1,0.24,0.04,1,0.0,42,0);

  //g_Stat1->SetPoint(0,0,-1.0);
  //g_Stat2->SetPoint(0,0,-1.0);
  //g_Sys1->SetPoint(0,0,-1.0);
  //g_Sys2->SetPoint(0,0,-1.0);
  //g_Stat1->SetPoint(1,0,-1.0);
  //g_Stat2->SetPoint(1,0,-1.0);
  //g_Sys1->SetPoint(1,0,-1.0);
  //g_Sys2->SetPoint(1,0,-1.0);

  for(int i = 0; i <14; i++)
  {
    //double pt1st, pt1sy, rho1st, rho1sy;
    double pt2st, pt2sy, rho2st, rho2sy;
    g_Stat2->GetPoint(i,pt2st,rho2st);
    g_Sys2->GetPoint(i,pt2sy,rho2sy);
    g_Stat2->SetPoint(i,pt2st+0.03,rho2st); 
    g_Sys2->SetPoint(i,pt2sy+0.03,rho2sy); 
  }
  
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_Stat2,style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_Stat1,style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);

//  plotSysErrorsBox(g_Sys2,color_phi_2nd);//,energy);
//  plotSysErrorsBox(g_Sys1,color_phi_1st);//,energy);

  //g_Stat1->SetMarkerStyle(style_phi_1at);
  //g_Stat1->SetMarkerColor(color_phi_1st);
  //g_Stat1->SetLineColor(color_phi_1st);
  //g_Stat1->SetMarkerSize(1.3);

  //g_Sys1->SetMarkerStyle(20);
  //g_Sys1->SetMarkerSize(1.3);
  //g_Sys1->SetMarkerColor(kBlue);
  //g_Sys1->SetLineColor(kBlue);
  //g_Stat1->Draw("pE  same");

  //plotSysErrorsBox(g_Sys1,kBlue);//,energy);

  //g_Stat2->SetMarkerStyle(20);
  //g_Stat2->SetMarkerColor(kOrange+7);
  //g_Stat2->SetLineColor(kOrange+7);
  //g_Stat2->SetMarkerSize(1.3);

  //g_Sys2->SetMarkerStyle(20);
  //g_Sys2->SetMarkerSize(1.3);
  //g_Sys2->SetMarkerColor(kOrange+7);
  //g_Sys2->SetLineColor(kOrange+7);
  //g_Stat2->Draw("pE  same");
  //c_rho_SysError->SaveAs("figures/phi_ordercomp.eps");
  c_rho_SysError->SaveAs(Form("figures/phi_ordercomp_rapidity_%s.pdf",vmsa::mBeamEnergy[energy].c_str()));
  c_rho_SysError->SaveAs(Form("figures/phi_ordercomp_rapidity_%s.png",vmsa::mBeamEnergy[energy].c_str()));

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

    bSys[i_pt] = new TBox(pt-0.02,rho-err,pt+0.02,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}


