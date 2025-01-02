#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

using namespace std;

void calculateEffPID_t(const int energy = 4, const string date = "20240702")
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  string InPutFile = Form("data/file_PID_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File = TFile::Open(InPutFile.c_str());

  TH2D *h_mPID[10][2][20][20]; // x axis is nsigma_kaon, y axis is m^2
  TH1D *h_mPID_nsig[10][2][20][20];
  TH1D *h_mPID_m2[10][2][20][20]; 
  TF2 *f_mPID[10][2][20][20];
  TF1 *f_mPID_nsig[10][2][20][20];
  TF1 *f_mPID_m2[10][2][20][20];
  TF1 *f_1D_mPID[3][10][2][20][20]; // first array runs over particle type: pion, kaon, proton
  TF1 *f_1D_mPID_nsig[3][10][2][20][20];
  TF1 *f_1D_mPID_m2[3][10][2][20][20];
  TGraphAsymmErrors *g_mPID_m2_pt_norm[3][10][2][20];
  TGraphAsymmErrors *g_mPID_m2_pt_nu[3][10][2][20];
  TGraphAsymmErrors *g_mPID_m2_pt_mean[3][10][2][20];
  TGraphAsymmErrors *g_mPID_m2_pt_sigma[3][10][2][20];

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 1; ic < 2; ic++)
    {
      for(int ip = 0; ip < 3; ip++)
      {
        for(int ieta = 0; ieta < 20; ieta++)
        {
          g_mPID_m2_pt_norm[ip][icent][ic][ieta]  = new TGraphAsymmErrors();
          g_mPID_m2_pt_nu[ip][icent][ic][ieta]  = new TGraphAsymmErrors();
          g_mPID_m2_pt_mean[ip][icent][ic][ieta]  = new TGraphAsymmErrors();
          g_mPID_m2_pt_sigma[ip][icent][ic][ieta] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = 0; ipt < 20; ipt++)
      {
        for(int ieta = 10; ieta < 11; ieta++)
        {
          cout << "pt = " << ipt << endl;
          string KEY_InPut = Form("h_mPID_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
          h_mPID[icent][ic][ipt][ieta] = (TH2D*)File->Get(KEY_InPut.c_str())->Clone();  
 
          string KEY_InPut_nsig = Form("h_mPID_nsig_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
          h_mPID_nsig[icent][ic][ipt][ieta] = h_mPID[icent][ic][ipt][ieta]->ProjectionX(KEY_InPut_nsig.c_str(),0,-1,"e");
 
          string KEY_InPut_m2 = Form("h_mPID_m2_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
          h_mPID_m2[icent][ic][ipt][ieta] = h_mPID[icent][ic][ipt][ieta]->ProjectionY(KEY_InPut_m2.c_str(),0,-1,"e");

          string FuncName = Form("f_mPID_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
          //f_mPID[icent][ic][ipt][ieta] = new TF2(FuncName.c_str(),ThreeParticleGaussian2D,-10.0,10.0,-0.5,1.5,18);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(0,h_mPID[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(1,0.0);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(2,0.2);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(3,h_mPID[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(4,0.0);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(5,0.2);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(6,h_mPID[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(7,0.0);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(8,0.2);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(9,h_mPID[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(10,0.01822044);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(11,0.2);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(12,h_mPID[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(13,0.247618);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(14,0.2);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(15,h_mPID[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(16,0.880354);
          //f_mPID[icent][ic][ipt][ieta]->SetParameter(17,0.2);

          //h_mPID[icent][ic][ipt][ieta]->Fit(f_mPID[icent][ic][ipt][ieta],"NMRI");

          FuncName = Form("f_mPID_m2_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
          f_mPID_m2[icent][ic][ipt][ieta] = new TF1(FuncName.c_str(),generalized_t_pdf_3part,-0.5,1.5,12);
          f_mPID_m2[icent][ic][ipt][ieta]->SetNpx(1000);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum());
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(1,100);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(2,0.01822044);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParLimits(2,0.01822044-0.1,0.01822044+0.1);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(3,0.002+0.006*ipt);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(4,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum()/3.0);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(5,100);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(6,0.247618);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParLimits(6,0.247618-0.1,0.247618+0.1);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(7,0.008+0.012*ipt);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(8,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum()/2.0);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(9,100);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(10,0.880354);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParLimits(10,0.880354-0.2,0.880354+0.2);
          f_mPID_m2[icent][ic][ipt][ieta]->SetParameter(11,0.044+0.012*ipt);

          h_mPID_m2[icent][ic][ipt][ieta]->Fit(f_mPID_m2[icent][ic][ipt][ieta],"NMR");


          for(int ip = 0; ip < 3; ip++)
          {
            g_mPID_m2_pt_norm[ip][icent][ic][ieta]->SetPoint(ipt,0.1+ipt*0.2,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(4*ip));
            g_mPID_m2_pt_norm[ip][icent][ic][ieta]->SetPointError(ipt,0.0,0.0,f_mPID_m2[icent][ic][ipt][ieta]->GetParError(4*ip),f_mPID_m2[icent][ic][ipt][ieta]->GetParError(4*ip));

            g_mPID_m2_pt_nu[ip][icent][ic][ieta]->SetPoint(ipt,0.1+ipt*0.2,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(1+4*ip));
            g_mPID_m2_pt_nu[ip][icent][ic][ieta]->SetPointError(ipt,0.0,0.0,f_mPID_m2[icent][ic][ipt][ieta]->GetParError(1+4*ip),f_mPID_m2[icent][ic][ipt][ieta]->GetParError(1+4*ip));

            g_mPID_m2_pt_mean[ip][icent][ic][ieta]->SetPoint(ipt,0.1+ipt*0.2,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(2+4*ip));
            g_mPID_m2_pt_mean[ip][icent][ic][ieta]->SetPointError(ipt,0.0,0.0,f_mPID_m2[icent][ic][ipt][ieta]->GetParError(2+4*ip),f_mPID_m2[icent][ic][ipt][ieta]->GetParError(2+4*ip));

            g_mPID_m2_pt_sigma[ip][icent][ic][ieta]->SetPoint(ipt,0.1+ipt*0.2,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(3+4*ip));
            g_mPID_m2_pt_sigma[ip][icent][ic][ieta]->SetPointError(ipt,0.0,0.0,f_mPID_m2[icent][ic][ipt][ieta]->GetParError(3+4*ip),f_mPID_m2[icent][ic][ipt][ieta]->GetParError(3+4*ip));
          }
          //FuncName = Form("f_mPID_nsig_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
          //f_mPID_nsig[icent][ic][ipt][ieta] = new TF1(FuncName.c_str(),ThreeParticleGaussian,-10.0,10.0,9);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(0,h_mPID_nsig[icent][ic][ipt][ieta]->GetMaximum());
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(1,0.01822044);
          ////f_mPID_nsig[icent][ic][ipt][ieta]->SetParLimits(1,0.01822044-0.05,0.01822044+0.05);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(2,0.1);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(3,h_mPID_nsig[icent][ic][ipt][ieta]->GetMaximum()/3.0);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(4,0.247618);
          ////f_mPID_nsig[icent][ic][ipt][ieta]->SetParLimits(4,0.247618-0.05,0.247618+0.05);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(5,0.1);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(6,h_mPID_nsig[icent][ic][ipt][ieta]->GetMaximum()/2.0);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(7,0.880354);
          ////f_mPID_nsig[icent][ic][ipt][ieta]->SetParLimits(7,0.880354-0.05,0.880354+0.05);
          //f_mPID_nsig[icent][ic][ipt][ieta]->SetParameter(8,0.1);

          //h_mPID_nsig[icent][ic][ipt][ieta]->Fit(f_mPID_nsig[icent][ic][ipt][ieta],"NMRI");
        }
      }
    }
  }

  string particlename[3] = {"pion","kaon","proton"};
  string charge[2] = {"K-","K+"};
  int cent_low[10]  = {70,60,50,40,30,20,10,5,0,20};
  int cent_high[10] = {80,70,60,50,40,30,20,10,5,60};

  TCanvas *cparam = new TCanvas("c","c",10,10,2000,1600);
  cparam->Divide(5,4);
  for(int i = 0; i < 20; i++)
  {
    cparam->cd(i+1);
    cparam->cd(i+1)->SetLeftMargin(0.15);
    cparam->cd(i+1)->SetBottomMargin(0.15);
    cparam->cd(i+1)->SetTicks(1,1);
    cparam->cd(i+1)->SetGrid(0,0);
  }

  string outputname = Form("figures/%s/EfficiencyPID_m2_1D_means_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop  = Form("%s]",outputname.c_str());

  cparam->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 1; ic < 2; ic++)
    {
      for(int ipart = 0; ipart < 3; ipart++)
      {
        for(int ieta = 10; ieta < 11; ieta++)
        {
          cparam->cd(ieta+1);
          g_mPID_m2_pt_mean[ipart][icent][ic][ieta]->SetStats(0);
          //double pt_low  = 0.1 + 0.2*double(ipt);
          //double pt_high = 0.1 + 0.2*double(ipt+1);
          g_mPID_m2_pt_mean[ipart][icent][ic][ieta]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-10)/10.,float(ieta-9)/10.));
          g_mPID_m2_pt_mean[ipart][icent][ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          g_mPID_m2_pt_mean[ipart][icent][ic][ieta]->GetYaxis()->SetTitle(Form("#mu_{%s}",particlename[ipart].c_str()));
          g_mPID_m2_pt_mean[ipart][icent][ic][ieta]->SetMarkerStyle(24);
          g_mPID_m2_pt_mean[ipart][icent][ic][ieta]->Draw("APE");

          cparam->Update();
          cparam->Print(outputname.c_str());
        }
      }
    }
  }
  cparam->Print(outputstop.c_str());

  outputname = Form("figures/%s/EfficiencyPID_m2_1D_sigmas_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 1; ic < 2; ic++)
    {
      for(int ipart = 0; ipart < 3; ipart++)
      {
        for(int ieta = 10; ieta < 11; ieta++)
        {
          cparam->cd(ieta+1);
          g_mPID_m2_pt_sigma[ipart][icent][ic][ieta]->SetStats(0);
          //double pt_low  = 0.1 + 0.2*double(ipt);
          //double pt_high = 0.1 + 0.2*double(ipt+1);
          g_mPID_m2_pt_sigma[ipart][icent][ic][ieta]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-10)/10.,float(ieta-9)/10.));
          g_mPID_m2_pt_sigma[ipart][icent][ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          g_mPID_m2_pt_sigma[ipart][icent][ic][ieta]->GetYaxis()->SetTitle(Form("#sigma_{%s}",particlename[ipart].c_str()));
          g_mPID_m2_pt_sigma[ipart][icent][ic][ieta]->SetMarkerStyle(24);
          g_mPID_m2_pt_sigma[ipart][icent][ic][ieta]->Draw("APE");

          cparam->Update();
          cparam->Print(outputname.c_str());
        }
      }
    }
  }
  cparam->Print(outputstop.c_str());

  outputname = Form("figures/%s/EfficiencyPID_m2_1D_nus_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 1; ic < 2; ic++)
    {
      for(int ipart = 0; ipart < 3; ipart++)
      {
        for(int ieta = 10; ieta < 11; ieta++)
        {
          cparam->cd(ieta+1);
          g_mPID_m2_pt_nu[ipart][icent][ic][ieta]->SetStats(0);
          //double pt_low  = 0.1 + 0.2*double(ipt);
          //double pt_high = 0.1 + 0.2*double(ipt+1);
          g_mPID_m2_pt_nu[ipart][icent][ic][ieta]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-10)/10.,float(ieta-9)/10.));
          g_mPID_m2_pt_nu[ipart][icent][ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          g_mPID_m2_pt_nu[ipart][icent][ic][ieta]->GetYaxis()->SetTitle(Form("#nu_{%s}",particlename[ipart].c_str()));
          g_mPID_m2_pt_nu[ipart][icent][ic][ieta]->SetMarkerStyle(24);
          g_mPID_m2_pt_nu[ipart][icent][ic][ieta]->Draw("APE");

          cparam->Update();
          cparam->Print(outputname.c_str());
        }
      }
    }
  }
  cparam->Print(outputstop.c_str());

  outputname = Form("figures/%s/EfficiencyPID_m2_1D_norms_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 1; ic < 2; ic++)
    {
      for(int ipart = 0; ipart < 3; ipart++)
      {
        for(int ieta = 10; ieta < 11; ieta++)
        {
          cparam->cd(ieta+1);
          g_mPID_m2_pt_norm[ipart][icent][ic][ieta]->SetStats(0);
          //double pt_low  = 0.1 + 0.2*double(ipt);
          //double pt_high = 0.1 + 0.2*double(ipt+1);
          g_mPID_m2_pt_norm[ipart][icent][ic][ieta]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-10)/10.,float(ieta-9)/10.));
          g_mPID_m2_pt_norm[ipart][icent][ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          g_mPID_m2_pt_norm[ipart][icent][ic][ieta]->GetYaxis()->SetTitle(Form("N_{%s}",particlename[ipart].c_str()));
          g_mPID_m2_pt_norm[ipart][icent][ic][ieta]->SetMarkerStyle(24);
          g_mPID_m2_pt_norm[ipart][icent][ic][ieta]->Draw("APE");

          cparam->Update();
          cparam->Print(outputname.c_str());
        }
      }
    }
  }
  cparam->Print(outputstop.c_str());


  TCanvas *c = new TCanvas("c","c",10,10,2000,1600);
  c->Divide(5,4);
  for(int i = 0; i < 20; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_nsigma_1D_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  c->Print(outputstart.c_str());

  //for(int icent = 0; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = 0; ipt < 20; ipt++)
  //    {
  //      for(int ieta = 0; ieta < 20; ieta++)
  //      {
  //        c->cd(ieta+1);
  //        h_mPID_nsig[icent][ic][ipt][ieta]->SetStats(0);
  //        double pt_low  = 0.1 + 0.2*double(ipt);
  //        double pt_high = 0.1 + 0.2*double(ipt+1);
  //        h_mPID_nsig[icent][ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high));
  //        h_mPID_nsig[icent][ic][ipt][ieta]->GetXaxis()->SetTitle("n#sigma_{K}");
  //        h_mPID_nsig[icent][ic][ipt][ieta]->GetYaxis()->SetTitle("Number of Tracks");
  //        h_mPID_nsig[icent][ic][ipt][ieta]->Draw("pE");
  //      }
  //      c->Update();
  //      c->Print(outputname.c_str());
  //    }
  //  }
  //}
  //c->Print(outputstop.c_str());


  outputname = Form("figures/%s/EfficiencyPID_m2_1D_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  int particlecolor[3] = {kOrange+7,kBlue,kGray+2};

  c->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 1; ic < 2; ic++)
    {
      for(int ipt = 0; ipt < 20; ipt++)
      {
        for(int ieta = 10; ieta < 11; ieta++)
        {
          c->cd(ieta+1);
          h_mPID_m2[icent][ic][ipt][ieta]->SetStats(0);
          double pt_low  = 0.1 + 0.2*double(ipt);
          double pt_high = 0.1 + 0.2*double(ipt+1);
          h_mPID_m2[icent][ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c, %1.1f<#eta<%1.1f",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-10)/10.,float(ieta-9)/10.));
          h_mPID_m2[icent][ic][ipt][ieta]->GetXaxis()->SetTitle("m^{2} (GeV/c^{2})");
          h_mPID_m2[icent][ic][ipt][ieta]->GetYaxis()->SetTitle("Number of Tracks");
          //h_mPID_m2[icent][ic][ipt][ieta]->SetMarkerStyle(24);
          h_mPID_m2[icent][ic][ipt][ieta]->Draw("EX0");

          f_mPID_m2[icent][ic][ipt][ieta]->SetLineColor(kMagenta);
          f_mPID_m2[icent][ic][ipt][ieta]->SetNpx(1000);
          f_mPID_m2[icent][ic][ipt][ieta]->SetLineWidth(1);
          f_mPID_m2[icent][ic][ipt][ieta]->Draw("same");
      
          for(int ipart = 0; ipart < 3; ipart++)
          {
            string FuncName = Form("f_mPID_m2_%s_cent%d_charge%d_pt%d_eta%d",particlename[ipart].c_str(),icent,ic,ipt,ieta);
            f_1D_mPID_m2[ipart][icent][ic][ipt][ieta] = new TF1(FuncName.c_str(),SingleParticleGaussian,-0.5,1.5,3);
            f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetNpx(1000);
            f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetParameter(0,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(0+ipart*3));
            f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetParameter(1,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(1+ipart*3));
            f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetParameter(2,f_mPID_m2[icent][ic][ipt][ieta]->GetParameter(2+ipart*3));
	    f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetLineColor(particlecolor[ipart]);
	    f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetLineWidth(1);
	    f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->SetLineStyle(1);
            f_1D_mPID_m2[ipart][icent][ic][ipt][ieta]->Draw("same");

            PlotLine(0.01822044-0.1,0.01822044-0.1,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum(),kBlack,1,2);
            PlotLine(0.01822044+0.1,0.01822044+0.1,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum(),kBlack,1,2);
            PlotLine(0.247618-0.1,0.247618-0.1,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum()    ,kBlack,1,2);
            PlotLine(0.247618+0.1,0.247618+0.1,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum()    ,kBlack,1,2);
            PlotLine(0.880354-0.2,0.880354-0.2,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum()    ,kBlack,1,2);
            PlotLine(0.880354+0.2,0.880354+0.2,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum()    ,kBlack,1,2);
 
            PlotLine(0.16,0.16,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum(),kGreen,1,2);
            PlotLine(0.36,0.36,0,h_mPID_m2[icent][ic][ipt][ieta]->GetMaximum(),kGreen,1,2);
          }
          
        }
        c->Update();
        c->Print(outputname.c_str());
      }
    }
  }
  c->Print(outputstop.c_str());


  //outputname = Form("figures/%s/EfficiencyPID_2D_studentt.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //c->Print(outputstart.c_str());

  //for(int icent = 0; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = 0; ipt < 20; ipt++)
  //    {
  //      for(int ieta = 0; ieta < 20; ieta++)
  //      {
  //        c->cd(ieta+1);
  //        c->cd(ieta+1)->SetRightMargin(0.15);
  //        h_mPID[icent][ic][ipt][ieta]->SetStats(0);
  //        double pt_low  = 0.1 + 0.2*double(ipt);
  //        double pt_high = 0.1 + 0.2*double(ipt+1);
  //        h_mPID[icent][ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high));
  //        h_mPID[icent][ic][ipt][ieta]->GetXaxis()->SetTitle("n#sigma_{K}");
  //        h_mPID[icent][ic][ipt][ieta]->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})");
  //        h_mPID[icent][ic][ipt][ieta]->GetZaxis()->SetTitle("Number of Tracks");
  //        h_mPID[icent][ic][ipt][ieta]->Draw("colz");
  //      }
  //      c->Update();
  //      c->Print(outputname.c_str());
  //    }
  //  }
  //}
  //c->Print(outputstop.c_str());
}
