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

void compRhoEff_CentEtaStudy(Int_t energy = 3, Int_t pid = 0, std::string setting = "Raw", int order = 1)
{
  std::string ep = "";
  if(order == 1) ep = "_FirstOrder";
  std::string EP[2] = {"","2nd"};
  string inputfileF0 = Form("../output/AuAu%s/%s/RhoCent_%sSysErrors_eta0p4%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),setting.c_str(),ep.c_str());
  string inputfileF1 = Form("../output/AuAu%s/%s/RhoCent_%sSysErrors_eta0p6%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),setting.c_str(),ep.c_str());
  string inputfileF2 = Form("../output/AuAu%s/%s/RhoCent_%sSysErrors_eta0p8%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),setting.c_str(),ep.c_str());
  string inputfileF3 = Form("../output/AuAu%s/%s/RhoCent_%sSysErrors_eta1_eta1%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),setting.c_str(),ep.c_str());

  TFile *File_InPutF0 = TFile::Open(inputfileF0.c_str());
  TFile *File_InPutF1 = TFile::Open(inputfileF1.c_str());
  TFile *File_InPutF2 = TFile::Open(inputfileF2.c_str());
  TFile *File_InPutF3 = TFile::Open(inputfileF3.c_str());

  TString StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TString SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TGraphAsymmErrors *g_StatF0   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF0->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_SysF0    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF0->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_StatF1   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF1->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_SysF1    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF1->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_StatF2   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF2->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_SysF2    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF2->Get( SysErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_StatF3   = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF3->Get(StatErrorRho.Data()))->Clone();
  TGraphAsymmErrors *g_SysF3    = (TGraphAsymmErrors*)((TGraphAsymmErrors*)   File_InPutF3->Get( SysErrorRho.Data()))->Clone();
 
  TGraMap g_mRhoF0_Stat;
  TGraMap g_mRhoF1_Stat;
  TGraMap g_mRhoF2_Stat;
  TGraMap g_mRhoF3_Stat;
  TGraMap g_mRhoF0_Sys;
  TGraMap g_mRhoF1_Sys;
  TGraMap g_mRhoF2_Sys;
  TGraMap g_mRhoF3_Sys;
  for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
  {           
    string KEY_rho = Form("rhoRaw_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    string KEY_rhoStat = Form("rhoRawStat_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    string KEY_rhoSys  = Form("rhoRawSys_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    g_mRhoF0_Stat[KEY_rho] = (TGraphAsymmErrors*)File_InPutF0->Get(KEY_rhoStat.c_str())->Clone();
    g_mRhoF1_Stat[KEY_rho] = (TGraphAsymmErrors*)File_InPutF1->Get(KEY_rhoStat.c_str())->Clone();
    g_mRhoF2_Stat[KEY_rho] = (TGraphAsymmErrors*)File_InPutF2->Get(KEY_rhoStat.c_str())->Clone();
    g_mRhoF3_Stat[KEY_rho] = (TGraphAsymmErrors*)File_InPutF3->Get(KEY_rhoStat.c_str())->Clone();
    g_mRhoF0_Sys[KEY_rho]  = (TGraphAsymmErrors*)File_InPutF0->Get(KEY_rhoSys.c_str())->Clone();
    g_mRhoF1_Sys[KEY_rho]  = (TGraphAsymmErrors*)File_InPutF1->Get(KEY_rhoSys.c_str())->Clone();
    g_mRhoF2_Sys[KEY_rho]  = (TGraphAsymmErrors*)File_InPutF2->Get(KEY_rhoSys.c_str())->Clone();
    g_mRhoF3_Sys[KEY_rho]  = (TGraphAsymmErrors*)File_InPutF3->Get(KEY_rhoSys.c_str())->Clone();
  }

  /*for(int i = 0; i < g_StatF0->GetN(); i++)
  {
    double pt, rho;
    g_StatF0->GetPoint(i,pt,rho);
    g_StatF0->SetPoint(i,pt-0.05,rho);
    g_StatF1->GetPoint(i,pt,rho);
    g_StatF1->SetPoint(i,pt+0.05,rho);
    //g_StatF2->GetPoint(i,pt,rho);
    //g_StatF2->SetPoint(i,pt+0.27,rho);

    g_SysF0->GetPoint(i,pt,rho);
    g_SysF0->SetPoint(i,pt-0.05,rho);
    g_SysF1->GetPoint(i,pt,rho);
    g_SysF1->SetPoint(i,pt+0.05,rho);
    //g_SysF2->GetPoint(i,pt,rho);
    //g_SysF2->SetPoint(i,pt+0.27,rho);

  }*/

  TH1F *h_frame = (TH1F*)File_InPutF1->Get("h_frame");

  TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //cout << "All good" << endl;

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->GetYaxis()->SetRangeUser(0.15,0.55);
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

  PlotLine(0.0,80.0,1.0/3.0,1.0/3.0,1,2,2);

  gStyle->SetErrorX(0);

  //besi19->SetLineColor(kBlack);
  //besi19->SetMarkerStyle(20);
  //besi19->SetMarkerSize(1.3);
  //besi19->SetMarkerColor(kBlack);
  //besi19->SetLineColor(kBlack);

  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->SetMarkerStyle(20);
  //besi19_sys->SetMarkerSize(1.3);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //plotSysErrorsBox(besi19_sys,kBlack);//,energy);
  //besi19->Draw("pE same");
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

  Draw_TGAE_Point_new_Symbol(5.0,0.48,0.0,0.0,0.0,0.0,20,kBlue,1.3); 
  string leg_count0 = Form("#phi |#eta|<0.4");
  plotTopLegend((char*)leg_count0.c_str(),6.0,0.477,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(5.0,0.46,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  string leg_count1 = Form("#phi |#eta|<0.6");
  plotTopLegend((char*)leg_count1.c_str(),6.0,0.457,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(5.0,0.44,0.0,0.0,0.0,0.0,20,kGray+2,1.3); 
  string leg_count2 = Form("#phi |#eta|<0.8");
  plotTopLegend((char*)leg_count2.c_str(),6.0,0.437,0.03,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(5.0,0.42,0.0,0.0,0.0,0.0,20,kBlack,1.3); 
  string leg_count3 = Form("#phi |#eta|<1.0");
  plotTopLegend((char*)leg_count3.c_str(),6.0,0.417,0.03,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,kOrange+7,1.3); 
  //string leg_countA = "Eff+Acc+Res Corrected #phi BES-II";
  //plotTopLegend((char*)leg_countA.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  //string leg_besi = "#phi BES-I";
  //Draw_TGAE_Point_new_Symbol(0.5,0.41,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  //plotTopLegend((char*)leg_besi.c_str(),0.6,0.407,0.03,1,0.0,42,0);

  //PlotLine(1,5,0.262,0.262,1,2,2);
  //string leg_line = "#rho_{00} = 1/3";
  //plotTopLegend((char*)leg_line.c_str(),0.85,0.26,0.03,1,0.0,42,0);


  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),30.0,0.28,0.04,1,0.0,42,0);
  //plotTopLegend((char*)"20%-60%",3.1,0.24,0.04,1,0.0,42,0);

  //g_StatF0->SetPoint(0,0,-1.0);
  //g_StatF1->SetPoint(0,0,-1.0);
  //g_StatF2->SetPoint(0,0,-1.0);
  //g_SysF0->SetPoint(0,0,-1.0);
  //g_SysF1->SetPoint(0,0,-1.0);
  //g_SysF2->SetPoint(0,0,-1.0);
  //g_StatF0->SetPoint(1,0,-1.0);
  //g_StatF1->SetPoint(1,0,-1.0);
  //g_StatF2->SetPoint(1,0,-1.0);
  //g_SysF0->SetPoint(1,0,-1.0);
  //g_SysF1->SetPoint(1,0,-1.0);
  //g_SysF2->SetPoint(1,0,-1.0);
  
  g_StatF0->SetMarkerStyle(20);
  g_StatF0->SetMarkerColor(kBlue);
  g_StatF0->SetLineColor(kBlue);
  g_StatF0->SetMarkerSize(1.3);

  g_SysF0->SetMarkerStyle(20);
  g_SysF0->SetMarkerSize(1.3);
  g_SysF0->SetMarkerColor(kBlue);
  g_SysF0->SetLineColor(kBlue);
  g_StatF0->Draw("pE  same");

  plotSysErrorsBox(g_SysF0,kBlue);//,energy);

  g_StatF1->SetMarkerStyle(20);
  g_StatF1->SetMarkerColor(kOrange+7);
  g_StatF1->SetLineColor(kOrange+7);
  g_StatF1->SetMarkerSize(1.3);

  g_SysF1->SetMarkerStyle(20);
  g_SysF1->SetMarkerSize(1.3);
  g_SysF1->SetMarkerColor(kOrange+7);
  g_SysF1->SetLineColor(kOrange+7);
  g_StatF1->Draw("pE  same");

  plotSysErrorsBox(g_SysF1,kOrange+7);//,energy);

  g_StatF2->SetMarkerStyle(20);
  g_StatF2->SetMarkerColor(kGray+2);
  g_StatF2->SetLineColor(kGray+2);
  g_StatF2->SetMarkerSize(1.3);

  g_SysF2->SetMarkerStyle(20);
  g_SysF2->SetMarkerSize(1.3);
  g_SysF2->SetMarkerColor(kGray+2);
  g_SysF2->SetLineColor(kGray+2);
  g_StatF2->Draw("pE  same");

  plotSysErrorsBox(g_SysF2,kGray+2);//,energy);

  g_StatF3->SetMarkerStyle(20);
  g_StatF3->SetMarkerColor(kBlack);
  g_StatF3->SetLineColor(kBlack);
  g_StatF3->SetMarkerSize(1.3);

  g_SysF3->SetMarkerStyle(20);
  g_SysF3->SetMarkerSize(1.3);
  g_SysF3->SetMarkerColor(kBlack);
  g_SysF3->SetLineColor(kBlack);
  g_StatF3->Draw("pE  same");

  plotSysErrorsBox(g_SysF3,kBlack);//,energy);

  string outputname = Form("./figures/Phi/19GeV/rho00Cent_allEta.pdf");
  string output_start = Form("%s[",outputname.c_str());

  TCanvas *c_pt = new TCanvas("c_pt","c_pt",10,10,800,800);
  c_pt->Divide(2,2);

  c_pt->Print(output_start.c_str());

  std::string centStrings[3] = {"40-80","10-40","0-10"};


  for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
  {
    c_pt->cd(i_pt+1);
    c_pt->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pt->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pt->cd(i_pt+1)->SetTicks(1,1);
    c_pt->cd(i_pt+1)->SetGrid(0,0);
    h_frame->SetTitle(Form("%.2f<p_{T}<%.2f GeV/c",vmsa::pt_low_cent[energy][i_pt],vmsa::pt_up_cent[energy][i_pt]));
    h_frame->DrawCopy("pE");
    PlotLine(0.0,80.0,1.0/3.0,1.0/3.0,1,2,2);
    //Draw_TGAE_Point_new_Symbol(-0.1,0.46,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 

    //plotTopLegend((char*)leg_count.c_str(),0.0,0.457,0.03,1,0.0,42,0);

    //PlotLine(-0.25,-0.05,0.438,0.438,1,2,2);
    //plotTopLegend((char*)leg_line.c_str(),0.0,0.437,0.03,1,0.0,42,0);

    //plotTopLegend((char*)leg_energy.c_str(),-0.3,0.305,0.04,1,0.0,42,0);
    //plotTopLegend((char*)centStrings[i_cent].c_str(),-0.1,0.295,0.04,1,0.0,42,0);
    //plotTopLegend((char*)Form("%.2f<p_{T}<%.2f GeV/c",vmsa::pt_low_y[energy][i_pt],vmsa::pt_up_y[energy][i_pt]),-0.3,0.285,0.04,1,0.0,42,0);


    string KEY_Default = Form("rhoRaw_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF0_Stat[KEY_Default],20,kBlue,0.0,1.0);
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF1_Stat[KEY_Default],20,kOrange+7,0.0,1.0);
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF2_Stat[KEY_Default],20,kGray+2,0.0,1.0);
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRhoF3_Stat[KEY_Default],20,kBlack,0.0,1.0);

    plotSysErrorsBox(g_mRhoF0_Sys[KEY_Default],kBlue);     
    plotSysErrorsBox(g_mRhoF1_Sys[KEY_Default],kOrange+7);     
    plotSysErrorsBox(g_mRhoF2_Sys[KEY_Default],kGray+2);     
    plotSysErrorsBox(g_mRhoF3_Sys[KEY_Default],kBlack);     

    
  }
  c_pt->Update();
  c_pt->Print(outputname.c_str());
  string output_stop = Form("%s]",outputname.c_str());
  c_pt->Print(output_stop.c_str()); // close pdf file


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

    bSys[i_pt] = new TBox(pt-1.6,rho-err,pt+1.6,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}
