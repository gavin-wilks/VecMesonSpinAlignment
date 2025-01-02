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
#include "TFitResultPtr.h"


using namespace std;

void calculateEffPID_phi_nsigtest(const int energy = 4, const string date = "20240703")
{
  int pt_max = 12;

  double kaonlimitlow = 0.075;
  double kaonlimitup = 0.15;

  const int neta = 10;
  const int nphi = 12;

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  string InPutFile = Form("data/file_PID_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile = Form("data/Phi_PIDEff_MoreEtaBins.root");
  TFile *File = TFile::Open(InPutFile.c_str());

  TH2D *h_mPID[2][20][neta][nphi]; // x axis is nsigma_kaon, y axis is m^2
  TH1D *h_mPID_nsig[2][20][neta][nphi];
  TH1D *h_mPID_m2[2][20][neta][nphi]; 
  TF2 *f_mPID[2][20][neta][nphi];
  TF1 *f_mPID_nsig[2][20][neta][nphi];
  TF1 *f_mPID_m2[2][20][neta][nphi];
  TFitResultPtr fitResult_mPID_m2[2][20][neta][nphi];
  TFitResultPtr fitResult_mPID_nsig[2][20][neta][nphi];
  TF1 *f_1D_mPID[3][2][20][neta][nphi]; // first array runs over particle type: pion, kaon, proton
  TF1 *f_1D_mPID_nsig[3][2][20][neta][nphi];
  TF1 *f_1D_mPID_m2[3][2][20][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_norm[3][2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_mean[3][2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_sigma[3][2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_ratio[2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_eta_ratio[2][20][nphi];
  TGraphAsymmErrors *g_mPID_m2_phi_ratio[2][20][neta];
  TGraphAsymmErrors *g_mPID_nsig_pt_ratio[2][neta][nphi];
  TGraphAsymmErrors *g_mPID_nsig_eta_ratio[2][20][nphi];
  TGraphAsymmErrors *g_mPID_nsig_phi_ratio[2][20][neta];

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ip = 0; ip < 3; ip++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          for(int iphi = 0; iphi < nphi; iphi++)
          {
            g_mPID_m2_pt_norm[ip][ic][ieta][iphi]  = new TGraphAsymmErrors();
            g_mPID_m2_pt_mean[ip][ic][ieta][iphi]  = new TGraphAsymmErrors();
            g_mPID_m2_pt_sigma[ip][ic][ieta][iphi] = new TGraphAsymmErrors();
          }
        }
      }
      for(int ipt = 0; ipt < pt_max; ipt++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          for(int iphi = 0; iphi < nphi; iphi++)
          {
            cout << "pt = " << ipt << endl;
            
            if(ipt < pt_max-1)
            {
              string KEY_InPut = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
              h_mPID[ic][ipt][ieta][iphi] = (TH2D*)File->Get(KEY_InPut.c_str())->Clone();   
            }
            if(ipt == pt_max-1) 
            {            
              string KEY_InPut = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
              h_mPID[ic][ipt][ieta][iphi] = (TH2D*)File->Get(KEY_InPut.c_str())->Clone();   

              for(int ipt2 = pt_max; ipt2 < 20; ipt2++)
              {
                string KEY_InPut2 = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt2,ieta,iphi);
                h_mPID[ic][ipt][ieta][iphi]->Add((TH2D*)File->Get(KEY_InPut2.c_str())->Clone());   
              }
            } 
          
            string KEY_InPut_nsig = Form("h_mPID_nsig_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            int y_bin_start = h_mPID[ic][ipt][ieta][iphi]->GetYaxis()->FindBin(0.16);
            int y_bin_stop  = h_mPID[ic][ipt][ieta][iphi]->GetYaxis()->FindBin(0.36);
            h_mPID_nsig[ic][ipt][ieta][iphi] = h_mPID[ic][ipt][ieta][iphi]->ProjectionX(KEY_InPut_nsig.c_str(),y_bin_start,y_bin_stop,"e");
            //h_mPID_nsig[ic][ipt][ieta][iphi] = h_mPID[ic][ipt][ieta][iphi]->ProjectionX(KEY_InPut_nsig.c_str(),0,-1,"e");
 
            //string KEY_InPut_m2 = Form("h_mPID_m2_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            //h_mPID_m2[ic][ipt][ieta][iphi] = h_mPID[ic][ipt][ieta][iphi]->ProjectionY(KEY_InPut_m2.c_str(),0,-1,"e");

            //string FuncName = Form("f_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            //f_mPID[ic][ipt][ieta][iphi] = new TF2(FuncName.c_str(),ThreeParticleGaussian2D,-10.0,10.0,-0.5,1.5,18);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(0,h_mPID[ic][ipt][ieta][iphi]->GetMaximum());
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(1,0.0);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(2,0.2);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(3,h_mPID[ic][ipt][ieta][iphi]->GetMaximum());
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(4,0.0);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(5,0.2);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(6,h_mPID[ic][ipt][ieta][iphi]->GetMaximum());
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(7,0.0);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(8,0.2);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(9,h_mPID[ic][ipt][ieta][iphi]->GetMaximum());
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(10,0.01822044);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(11,0.2);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(12,h_mPID[ic][ipt][ieta][iphi]->GetMaximum());
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(13,0.247618);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(14,0.2);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(15,h_mPID[ic][ipt][ieta][iphi]->GetMaximum());
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(16,0.880354);
            //f_mPID[ic][ipt][ieta][iphi]->SetParameter(17,0.2);

            ////h_mPID[ic][ipt][ieta]->Fit(f_mPID[ic][ipt][ieta],"NMRI");

            string FuncName = Form("f_mPID_nsig_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            if(ipt < 5) 
            {
              f_mPID_nsig[ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),SingleParticleGaussian,-7,7,3);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(0,h_mPID_nsig[ic][ipt][ieta][iphi]->GetMaximum());
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParLimits(0,0.0,1000000000);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(1,0.0);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParLimits(1,-5,5);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(2,1.5);
            }
            else 
            {
              f_mPID_nsig[ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),TwoParticleGaussian,-7,7,6);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(0,h_mPID_nsig[ic][ipt][ieta][iphi]->GetMaximum()/1.5);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParLimits(0,0.0,1000000000);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(1,0.0);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParLimits(1,-1,0.5);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(2,1.5);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(3,h_mPID_nsig[ic][ipt][ieta][iphi]->GetMaximum()/3);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParLimits(3,0.0,1000000000);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(4,2.0);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParLimits(4,0.0,5);
              f_mPID_nsig[ic][ipt][ieta][iphi]->SetParameter(5,1.5);
            }

            //if(ipt < 2) 
            //{
            //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(6,0.0);
            //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(7,0.0);
            //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(8,0.0);
            //  if(ipt < 1) 
            //  {
            //    f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(3,0.0);
            //    f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(4,0.0);
            //    f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(5,0.0);
            //  }
            //}

            //double ptmean = double(ipt)*0.2;
            //if(ipt == 13) ptmean = 3.4;

            //if(ipt > 8)  
            //{
            //  TF1 *tempN[3];           
            //  TF1 *tempSigma[3];

            //  for(int ip = 0; ip < 3; ip++)
            //  {
            //    tempN[ip] = new TF1("tempN","[0]*exp(-[1]*x)",0.5,1.7);
            //    tempSigma[ip] = new TF1("tempSigma","[0]+[1]*x+[2]*x*x",0.5,1.7);

            //    double pt, val; 
            //    if(ip == 0) g_mPID_m2_pt_norm[0][ic][ieta][iphi]->GetPoint(0,pt,val);
            //    if(ip == 1) g_mPID_m2_pt_norm[1][ic][ieta][iphi]->GetPoint(1,pt,val);
            //    if(ip == 2) g_mPID_m2_pt_norm[2][ic][ieta][iphi]->GetPoint(2,pt,val);
            //    tempN[ip]->SetParameter(0,val);
 
            //    g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->Fit(tempN[ip],"NMRI");
            //    g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->Fit(tempSigma[ip],"NMRI");

            //  }
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(0,tempN[0]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(2,tempSigma[0]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(3,tempN[1]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(5,tempSigma[1]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(6,tempN[2]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(8,tempSigma[2]->Eval(ptmean));
            //}
          

            //fitResult_mPID_m2[ic][ipt][ieta][iphi] = h_mPID_m2[ic][ipt][ieta][iphi]->Fit(f_mPID_m2[ic][ipt][ieta][iphi],"NRIBS");


            //for(int ip = 0; ip < 3; ip++)
            //{

            //  g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->SetPoint(ipt,ptmean,  f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(3*ip));
            //  g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->SetPointError(ipt,0.0,0.0, f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(3*ip),  f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(3*ip));
            //  g_mPID_m2_pt_mean[ip][ic][ieta][iphi]->SetPoint(ipt,ptmean,  f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(1+3*ip));
            //  g_mPID_m2_pt_mean[ip][ic][ieta][iphi]->SetPointError(ipt,0.0,0.0, f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(1+3*ip),f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(1+3*ip));
            //  g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->SetPoint(ipt,ptmean, f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(2+3*ip));
            //  g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->SetPointError(ipt,0.0,0.0,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(2+3*ip),f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(2+3*ip));
            //}
            ////FuncName = Form("f_mPID_nsig_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
            ////f_mPID_nsig[ic][ipt][ieta] = new TF1(FuncName.c_str(),ThreeParticleGaussian,-10.0,10.0,9);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(0,h_mPID_nsig[ic][ipt][ieta]->GetMaximum());
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(1,0.01822044);
            //////f_mPID_nsig[ic][ipt][ieta]->SetParLimits(1,0.01822044-0.05,0.01822044+0.05);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(2,0.1);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(3,h_mPID_nsig[ic][ipt][ieta]->GetMaximum()/3.0);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(4,0.247618);
            //////f_mPID_nsig[ic][ipt][ieta]->SetParLimits(4,0.247618-0.05,0.247618+0.05);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(5,0.1);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(6,h_mPID_nsig[ic][ipt][ieta]->GetMaximum()/2.0);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(7,0.880354);
            //////f_mPID_nsig[ic][ipt][ieta]->SetParLimits(7,0.880354-0.05,0.880354+0.05);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParameter(8,0.1);

            fitResult_mPID_nsig[ic][ipt][ieta][iphi] = h_mPID_nsig[ic][ipt][ieta][iphi]->Fit(f_mPID_nsig[ic][ipt][ieta][iphi],"NRIBS");
          }
        }
      }
    }
  }

  string particlename[3] = {"pion","kaon","proton"};
  string charge[2] = {"K-","K+"};
  int cent_low[10]  = {70,60,50,40,30,20,10,5,0,20};
  int cent_high[10] = {80,70,60,50,40,30,20,10,5,60};

  TCanvas *cparam = new TCanvas("c","c",10,10,1600,1200);
  cparam->Divide(4,3);
  for(int i = 0; i < 12; i++)
  {
    cparam->cd(i+1);
    cparam->cd(i+1)->SetLeftMargin(0.15);
    cparam->cd(i+1)->SetBottomMargin(0.15);
    cparam->cd(i+1)->SetTicks(1,1);
    cparam->cd(i+1)->SetGrid(0,0);
  }

  string outputname = Form("figures/%s/EfficiencyPID_m2_1D_means.pdf",vmsa::mBeamEnergy[energy].c_str());
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop  = Form("%s]",outputname.c_str());

  //cparam->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipart = 0; ipart < 3; ipart++)
  //    {
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        for(int iphi = 0; iphi < nphi; iphi++)
  //        {
  //          cparam->cd(iphi+1);
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->SetStats(0);
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-5)/5.,float(ieta-4)/5.));
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->GetYaxis()->SetTitle(Form("#mu_{%s}",particlename[ipart].c_str()));
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->SetMarkerStyle(24);
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->Draw("APE");

  //        }
  //        cparam->Update();
  //        cparam->Print(outputname.c_str());
  //      }
  //    }
  //  }
  //}
  //cparam->Print(outputstop.c_str());

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_sigmas.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //cparam->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipart = 0; ipart < 3; ipart++)
  //    {
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        for(int iphi = 0; iphi < nphi; iphi++)
  //        {
  //          cparam->cd(iphi+1);
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->SetStats(0);
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-5)/5.,float(ieta-4)/5.));
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->GetYaxis()->SetTitle(Form("#sigma_{%s}",particlename[ipart].c_str()));
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->SetMarkerStyle(24);
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->Draw("APE");

  //        }
  //        cparam->Update();
  //        cparam->Print(outputname.c_str());
  //      }
  //    }
  //  }
  //}
  //cparam->Print(outputstop.c_str());

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_norms.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //cparam->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipart = 0; ipart < 3; ipart++)
  //    {
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        for(int iphi = 0; iphi < nphi; iphi++)
  //        {
  //          cparam->cd(iphi+1);
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->SetStats(0);
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-5)/5.,float(ieta-4)/5.));
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->GetYaxis()->SetTitle(Form("N_{%s}",particlename[ipart].c_str()));
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->SetMarkerStyle(24);
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->Draw("APE");

  //        }
  //        cparam->Update();
  //        cparam->Print(outputname.c_str());
  //      }
  //    }
  //  }
  //}
  //cparam->Print(outputstop.c_str());


  TCanvas *c = new TCanvas("c","c",10,10,1600,1200);
  c->Divide(4,3);
  for(int i = 0; i < 12; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_nsigma_1D.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  c->Print(outputstart.c_str());
  int particlecolor[3] = {kOrange+7,kBlue,kGray+2};

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          g_mPID_nsig_pt_ratio[ic][ieta][iphi] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = 0; ipt < pt_max; ipt++)
      {
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          g_mPID_nsig_eta_ratio[ic][ipt][iphi] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = 0; ipt < pt_max; ipt++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          g_mPID_nsig_phi_ratio[ic][ipt][ieta] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = 0; ipt < pt_max; ipt++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          for(int iphi = 0; iphi < nphi; iphi++)
          {
            c->cd(iphi+1);
            h_mPID_nsig[ic][ipt][ieta][iphi]->SetStats(0);
            double pt_low  = 0.1 + 0.2*double(ipt);
            double pt_high = 0.1 + 0.2*double(ipt+1);
            if(ipt == pt_max-1)
            {
              //pt_low = 2.7;
              pt_low = 4.1-(20-pt_max+1)*0.2;
              pt_high = 4.1;
            }
            double etamean = (double(ieta)-14.5)/15.0;
            double phimean = 0.0;//-TMath::Pi()*(double(iphi)-5.5)/6.0;

            h_mPID_nsig[ic][ipt][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c, %1.2f<#eta<%1.2f, %d#pi/6<#phi<%d#pi/6",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-neta/2)/(neta/2),float(ieta-neta/2)/(neta/2),iphi-6,iphi-5));
            h_mPID_nsig[ic][ipt][ieta][iphi]->GetXaxis()->SetTitle("n#sigma_{K}");
            h_mPID_nsig[ic][ipt][ieta][iphi]->GetYaxis()->SetTitle("Number of Tracks");
            h_mPID_nsig[ic][ipt][ieta][iphi]->Draw("pE");

            f_mPID_nsig[ic][ipt][ieta][iphi]->SetLineColor(kMagenta);
            f_mPID_nsig[ic][ipt][ieta][iphi]->SetNpx(1000);
            f_mPID_nsig[ic][ipt][ieta][iphi]->SetLineWidth(1);
            f_mPID_nsig[ic][ipt][ieta][iphi]->Draw("same");

            //int particlemax = 2;
            int particlemax = 1;
            if(ipt >= 5) particlemax = 2;
            for(int ipart = 0; ipart < particlemax; ipart++)
            {
              string FuncName = Form("f_mPID_nsig_%s_cent%d_charge%d_pt%d_eta%d_phi%d",particlename[ipart].c_str(),icent,ic,ipt,ieta,iphi);
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),SingleParticleGaussian,-7,7,3);
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetNpx(1000);
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetParameter(0,f_mPID_nsig[ic][ipt][ieta][iphi]->GetParameter(0+ipart*3));
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetParameter(1,f_mPID_nsig[ic][ipt][ieta][iphi]->GetParameter(1+ipart*3));
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetParameter(2,f_mPID_nsig[ic][ipt][ieta][iphi]->GetParameter(2+ipart*3));
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetParError(0,f_mPID_nsig[ic][ipt][ieta][iphi]->GetParError(0+ipart*3));
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetParError(1,f_mPID_nsig[ic][ipt][ieta][iphi]->GetParError(1+ipart*3));
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetParError(2,f_mPID_nsig[ic][ipt][ieta][iphi]->GetParError(2+ipart*3));

              if(ipart == 0)
              {
                double params[3] = {0.0};
                TMatrixDSym covMat(3);
                params[0] = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetParams()[0];
                params[1] = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetParams()[1];
                params[2] = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetParams()[2];
                covMat(0,0) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(0,0);
                covMat(0,1) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(0,1);
                covMat(0,2) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(0,2);
                covMat(1,0) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(1,0);
                covMat(1,1) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(1,1);
                covMat(1,2) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(1,2);
                covMat(2,0) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(2,0);
                covMat(2,1) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(2,1);
                covMat(2,2) = fitResult_mPID_nsig[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(2,2);

                float bin_width = h_mPID_nsig[ic][ipt][ieta][iphi]->GetBinWidth(1);

                double kaons_selected  = f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->Integral(-2.0,2.0)/bin_width;   
                double dkaons_selected = f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->IntegralError(-2.0,2.0,params,covMat.GetMatrixArray())/bin_width;   
                double kaons_total     = f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->Integral(-7,7)/bin_width;   
                double dkaons_total    = f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->IntegralError(-7,7,params,covMat.GetMatrixArray())/bin_width;
         
                if(kaons_total == 0.0) continue; 

                double covariance = dkaons_selected * dkaons_selected;

                cout << "kaons_selected = " << kaons_selected << " +/- " << dkaons_selected << endl;
                cout << "kaons_total    = " << kaons_total << " +/- " << dkaons_total << endl;
                
                double R = kaons_selected/kaons_total;
                if(R > 1+1e-6 || R < 0) continue;               

                double relError_selected = dkaons_selected/kaons_selected;   
                double relError_total    = dkaons_total   /kaons_total;   
                double correlatedTerm = -2.0*covariance/(kaons_selected*kaons_total);

                cout << "relError_selected = " << relError_selected << endl;
                cout << "relError_total    = " << relError_total << endl;
                cout << "correlatedTerm    = " << correlatedTerm << endl;
                double sqrtinnerds = relError_selected * relError_selected + relError_total * relError_total + correlatedTerm;
                cout << "sqrtinnerds = " << sqrtinnerds << endl;
                double deltaR = R * sqrt(sqrtinnerds);
                if(sqrtinnerds < 1e-10 && sqrtinnerds > -1e-10) deltaR = 0.0;
                if(deltaR > 1) continue;
               

                g_mPID_nsig_pt_ratio[ic][ieta][iphi]->SetPoint(ipt,(pt_low+pt_high)/2.0,R);
                g_mPID_nsig_pt_ratio[ic][ieta][iphi]->SetPointError(ipt,0.0,0.0,deltaR,deltaR);
                g_mPID_nsig_eta_ratio[ic][ipt][iphi]->SetPoint(ieta,etamean,R);
                g_mPID_nsig_eta_ratio[ic][ipt][iphi]->SetPointError(ieta,0.0,0.0,deltaR,deltaR);
                g_mPID_nsig_phi_ratio[ic][ipt][ieta]->SetPoint(iphi,phimean,R);
                g_mPID_nsig_phi_ratio[ic][ipt][ieta]->SetPointError(iphi,0.0,0.0,deltaR,deltaR);
                cout << "Filling plots with " << R << " +/- " << deltaR << endl;
              }

	      f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetLineColor(particlecolor[ipart]);
	      f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetLineWidth(1);
	      f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->SetLineStyle(1);
              f_1D_mPID_nsig[ipart][ic][ipt][ieta][iphi]->Draw("same");

              //PlotLine(0.01822044-0.1,0.01822044-0.1,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kBlack,1,2);
              //PlotLine(0.01822044+0.1,0.01822044+0.1,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kBlack,1,2);
              //PlotLine(0.247618-kaonlimitlow,0.247618-kaonlimitlow,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
              //PlotLine(0.247618+kaonlimitup,0.247618+kaonlimitup,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
              //PlotLine(0.880354-0.2,0.880354-0.2,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
              //PlotLine(0.880354+0.2,0.880354+0.2,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
 
              PlotLine(-2.0,-2.0,0,h_mPID_nsig[ic][ipt][ieta][iphi]->GetMaximum(),kGreen,1,2);
              PlotLine(2.0,2.0,0,h_mPID_nsig[ic][ipt][ieta][iphi]->GetMaximum(),kGreen,1,2);
            }
          }
          c->Update();
          c->Print(outputname.c_str());
        }
      }
    }
  }
  c->Print(outputstop.c_str());


  //outputname = Form("figures/%s/EfficiencyPID_m2_1D.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //int particlecolor[3] = {kOrange+7,kBlue,kGray+2};

  //c->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ieta = 0; ieta < neta; ieta++)
  //    {
  //      for(int iphi = 0; iphi < nphi; iphi++)
  //      {
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi] = new TGraphAsymmErrors();
  //      }
  //    }
  //    for(int ipt = 0; ipt < pt_max; ipt++)
  //    {
  //      for(int iphi = 0; iphi < nphi; iphi++)
  //      {
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi] = new TGraphAsymmErrors();
  //      }
  //    }
  //    for(int ipt = 0; ipt < 14; ipt++)
  //    {
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta] = new TGraphAsymmErrors();
  //      }
  //    }
  //    for(int ipt = 0; ipt < 14; ipt++)
  //    {
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        for(int iphi = 0; iphi < nphi; iphi++)
  //        {
  //          c->cd(iphi+1);
  //          h_mPID_m2[ic][ipt][ieta][iphi]->SetStats(0);
  //          double pt_low  = 0.1 + 0.2*double(ipt);
  //          double pt_high = 0.1 + 0.2*double(ipt+1);
  //          if(ipt == 13)
  //          {
  //            pt_low = 2.7;
  //            pt_high = 4.1;
  //          }
  //          double etamean = -(double(ieta)-4.5)/5.0;
  //          double phimean = -TMath::Pi()*(double(iphi)-5.5)/6.0;


  //          h_mPID_m2[ic][ipt][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c, %1.1f<#eta<%1.1f",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-5)/5.,float(ieta-4)/5.));
  //          h_mPID_m2[ic][ipt][ieta][iphi]->GetXaxis()->SetTitle("m^{2} (GeV/c^{2})");
  //          h_mPID_m2[ic][ipt][ieta][iphi]->GetYaxis()->SetTitle("Number of Tracks");
  //          h_mPID_m2[ic][ipt][ieta][iphi]->Draw("EX0");

  //          f_mPID_m2[ic][ipt][ieta][iphi]->SetLineColor(kMagenta);
  //          f_mPID_m2[ic][ipt][ieta][iphi]->SetNpx(1000);
  //          f_mPID_m2[ic][ipt][ieta][iphi]->SetLineWidth(1);
  //          f_mPID_m2[ic][ipt][ieta][iphi]->Draw("same");
  //    
  //          for(int ipart = 0; ipart < 3; ipart++)
  //          {
  //            string FuncName = Form("f_mPID_m2_%s_cent%d_charge%d_pt%d_eta%d_phi%d",particlename[ipart].c_str(),icent,ic,ipt,ieta,iphi);
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),SingleParticleGaussian,-0.5,1.5,3);
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetNpx(1000);
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParameter(0,f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(0+ipart*3));
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParameter(1,f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(1+ipart*3));
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParameter(2,f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(2+ipart*3));
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParError(0,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(0+ipart*3));
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParError(1,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(1+ipart*3));
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParError(2,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(2+ipart*3));

  //            if(ipart == 1 && ipt > 0)
  //            {
  //              double params[3] = {0.0};
  //              TMatrixDSym covMat(3);
  //              params[0] = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetParams()[3];
  //              params[1] = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetParams()[4];
  //              params[2] = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetParams()[5];
  //              covMat(0,0) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(3,3);
  //              covMat(0,1) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(3,4);
  //              covMat(0,2) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(3,5);
  //              covMat(1,0) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(4,3);
  //              covMat(1,1) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(4,4);
  //              covMat(1,2) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(4,5);
  //              covMat(2,0) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(5,3);
  //              covMat(2,1) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(5,4);
  //              covMat(2,2) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(5,5);

  //              float bin_width = h_mPID_m2[ic][ipt][ieta][iphi]->GetBinWidth(1);

  //              double kaons_selected  = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->Integral(0.16,0.36)/bin_width;   
  //              double dkaons_selected = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->IntegralError(0.16,0.36,params,covMat.GetMatrixArray())/bin_width;   
  //              double kaons_total     = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->Integral(-0.5,1.5)/bin_width;   
  //              double dkaons_total    = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->IntegralError(-0.5,1.5,params,covMat.GetMatrixArray())/bin_width;
  //       
  //              if(kaons_total == 0.0) continue; 

  //              double covariance = dkaons_selected * dkaons_selected;

  //              cout << "kaons_selected = " << kaons_selected << " +/- " << dkaons_selected << endl;
  //              cout << "kaons_total    = " << kaons_total << " +/- " << dkaons_total << endl;
  //              
  //              double R = kaons_selected/kaons_total;
  //              if(R > 1 || R < 0) continue;               

  //              double relError_selected = dkaons_selected/kaons_selected;   
  //              double relError_total    = dkaons_total   /kaons_total;   
  //              double correlatedTerm = -2.0*covariance/(kaons_selected*kaons_total);

  //              cout << "relError_selected = " << relError_selected << endl;
  //              cout << "relError_total    = " << relError_total << endl;
  //              cout << "correlatedTerm    = " << correlatedTerm << endl;
  //              double sqrtinnerds = relError_selected * relError_selected + relError_total * relError_total + correlatedTerm;
  //              cout << "sqrtinnerds = " << sqrtinnerds << endl;
  //              double deltaR = R * sqrt(sqrtinnerds);
  //              if(sqrtinnerds < 1e-10 && sqrtinnerds > -1e-10) deltaR = 0.0;
  //              if(deltaR > 1) continue;
  //             

  //              g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetPoint(ipt,(pt_low+pt_high)/2.0,R);
  //              g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetPointError(ipt,0.0,0.0,deltaR,deltaR);
  //              g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetPoint(ieta,etamean,R);
  //              g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetPointError(ieta,0.0,0.0,deltaR,deltaR);
  //              g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetPoint(iphi,phimean,R);
  //              g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetPointError(iphi,0.0,0.0,deltaR,deltaR);
  //              cout << "Filling plots with " << R << " +/- " << deltaR << endl;
  //            }

  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetLineColor(particlecolor[ipart]);
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetLineWidth(1);
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetLineStyle(1);
  //            f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->Draw("same");

  //            PlotLine(0.01822044-0.1,0.01822044-0.1,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kBlack,1,2);
  //            PlotLine(0.01822044+0.1,0.01822044+0.1,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kBlack,1,2);
  //            PlotLine(0.247618-kaonlimitlow,0.247618-kaonlimitlow,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
  //            PlotLine(0.247618+kaonlimitup,0.247618+kaonlimitup,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
  //            PlotLine(0.880354-0.2,0.880354-0.2,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
  //            PlotLine(0.880354+0.2,0.880354+0.2,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
 
  //            PlotLine(0.16,0.16,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kGreen,1,2);
  //            PlotLine(0.36,0.36,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kGreen,1,2);
  //          }
  //        } 
  //        c->Update();
  //        c->Print(outputname.c_str());
  //      }
  //    }
  //  }
  //}
  //c->Print(outputstop.c_str());

  //cout << "Finished with plots before ratios" << endl;

  //TCanvas *cparam2 = new TCanvas("c","c",10,10,1600,1200);
  //cparam2->Divide(4,3);
  //for(int i = 0; i < 12; i++)
  //{
  //  cparam2->cd(i+1);
  //  cparam2->cd(i+1)->SetLeftMargin(0.15);
  //  cparam2->cd(i+1)->SetBottomMargin(0.15);
  //  cparam2->cd(i+1)->SetTicks(1,1);
  //  cparam2->cd(i+1)->SetGrid(0,0);
  //}

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_pt_ratios.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //cparam2->Print(outputstart.c_str());

  //TF1 *m2_PID_pT[10][2][neta][nphi];

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ieta = 0; ieta < neta; ieta++)
  //    {
  //      for(int iphi = 0; iphi < nphi; iphi++)
  //      {
  //        cparam2->cd(iphi+1);
  //        //g_mPID_m2_pt_ratio[ic][ieta][iphi]->Print();
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetStats(0);
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<#eta<%1.1f, %d#pi/6<#phi<%d#pi/6",charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-5)/5.,float(ieta-4)/5.,iphi-6,iphi-5));
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->GetYaxis()->SetTitle("PID Efficiency");
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetMarkerStyle(24);
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->Draw("APE");
  //       
  //        //m2_PID_pT[icent][ic][ieta][iphi] = new TF1(Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi),m2_PID,0.3,2.5,5);
  //        m2_PID_pT[icent][ic][ieta][iphi] = new TF1(Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi),m2_PID_tanh,0.9,2.5,4);
  //        m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(0,0.25);
  //        m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(1,-5);
  //        m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(2,1.75);
  //        m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(3,0.75);
  //        //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(0,2);
  //        //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(1,3);
  //        //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(2,-7);
  //        //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(3,1.5);
  //        //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(4,0.5);
  //        g_mPID_m2_pt_ratio[ic][ieta][iphi]->Fit(m2_PID_pT[icent][ic][ieta][iphi],"NMRI");
  //        m2_PID_pT[icent][ic][ieta][iphi]->SetLineColor(kRed);
  //        m2_PID_pT[icent][ic][ieta][iphi]->Draw("same");
  //      }
  //      cparam2->Update();
  //      cparam2->Print(outputname.c_str());
  //    }
  //  }
  //}
  //cparam2->Print(outputstop.c_str());

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_eta_ratios.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //cparam2->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = 0; ipt < 14; ipt++)
  //    {
  //      double pt_low  = 0.1 + 0.2*double(ipt);
  //      double pt_high = 0.1 + 0.2*double(ipt+1);
  //      if(ipt == 13)
  //      {
  //        pt_low = 2.7;
  //        pt_high = 4.1;
  //      }
  //      for(int iphi = 0; iphi < nphi; iphi++)
  //      {
  //        cparam2->cd(iphi+1);
  //        //g_mPID_m2_eta_ratio[ic][ipt][iphi]->Print();
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetStats(0);
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<p_{T}<%1.1fGeV/c, %d#pi/6<#phi<%d#pi/6",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,iphi-6,iphi-5));
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi]->GetXaxis()->SetTitle("#eta");
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi]->GetYaxis()->SetTitle("PID Efficiency");
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetMarkerStyle(24);
  //        g_mPID_m2_eta_ratio[ic][ipt][iphi]->Draw("APE");
  //      }
  //      cparam2->Update();
  //      cparam2->Print(outputname.c_str());
  //    }
  //  }
  //}
  //cparam2->Print(outputstop.c_str());

  //TCanvas *cparam3 = new TCanvas("c","c",10,10,2000,800);
  //cparam3->Divide(5,2);
  //for(int i = 0; i < 10; i++)
  //{
  //  cparam3->cd(i+1);
  //  cparam3->cd(i+1)->SetLeftMargin(0.15);
  //  cparam3->cd(i+1)->SetBottomMargin(0.15);
  //  cparam3->cd(i+1)->SetTicks(1,1);
  //  cparam3->cd(i+1)->SetGrid(0,0);
  //}

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_phi_ratios.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //cparam3->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = 0; ipt < 14; ipt++)
  //    {
  //      double pt_low  = 0.1 + 0.2*double(ipt);
  //      double pt_high = 0.1 + 0.2*double(ipt+1);
  //      if(ipt == 13)
  //      {
  //        pt_low = 2.7;
  //        pt_high = 4.1;
  //      }
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        cparam3->cd(ieta+1);
  //        //g_mPID_m2_phi_ratio[ic][ipt][ieta]->Print();
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetStats(0);
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<p_{T}<%1.1fGeV/c, %1.1f<#eta<%1.1f",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-5)/5.,float(ieta-4)/5.));
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta]->GetXaxis()->SetTitle("#phi");
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta]->GetYaxis()->SetTitle("PID Efficiency");
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetMarkerStyle(24);
  //        g_mPID_m2_phi_ratio[ic][ipt][ieta]->Draw("APE");
  //      }
  //      cparam3->Update();
  //      cparam3->Print(outputname.c_str());
  //    }
  //  }
  //}
  //cparam3->Print(outputstop.c_str());

  //outputname = Form("figures/%s/EfficiencyPID_2D.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //c->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = 0; ipt < 14; ipt++)
  //    {
  //      for(int ieta = 0; ieta < neta; ieta++)
  //      {
  //        for(int iphi = 0; iphi < nphi; iphi++)
  //        {
  //          c->cd(iphi+1);
  //          c->cd(iphi+1)->SetRightMargin(0.15);
  //          h_mPID[ic][ipt][ieta][iphi]->SetStats(0);
  //          double pt_low  = 0.1 + 0.2*double(ipt);
  //          double pt_high = 0.1 + 0.2*double(ipt+1);
  //          h_mPID[ic][ipt][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c, %1.1f<#eta<%1.1f, %d#pi/6<#phi<%d#pi/6",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-5)/5.,float(ieta-4)/5.,iphi-6,iphi-5));
  //          h_mPID[ic][ipt][ieta][iphi]->GetXaxis()->SetTitle("n#sigma_{K}");
  //          h_mPID[ic][ipt][ieta][iphi]->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})");
  //          h_mPID[ic][ipt][ieta][iphi]->GetZaxis()->SetTitle("Number of Tracks");
  //          h_mPID[ic][ipt][ieta][iphi]->Draw("colz");
  //          PlotLine(-2.5,2.5,0.16,0.16, kRed,1,2);
  //          PlotLine(-2.5,2.5,0.36,0.36, kRed,1,2);
  //          PlotLine(-2.5,-2.5,0.16,0.36,kRed,1,2);
  //          PlotLine(2.5,2.5,0.16,0.36,  kRed,1,2);
  //        }
  //        c->Update();
  //        c->Print(outputname.c_str());
  //      }
  //    }
  //  }
  //}
  //c->Print(outputstop.c_str());

  TCanvas *cparam2 = new TCanvas("c","c",10,10,1600,1200);
  cparam2->Divide(4,3);
  for(int i = 0; i < 12; i++)
  {
    cparam2->cd(i+1);
    cparam2->cd(i+1)->SetLeftMargin(0.15);
    cparam2->cd(i+1)->SetBottomMargin(0.15);
    cparam2->cd(i+1)->SetTicks(1,1);
    cparam2->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_nsig_1D_pt_ratios.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam2->Print(outputstart.c_str());

  TF1 *nsig_PID_pT[10][2][neta][nphi];

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int iphi = 0; iphi < nphi; iphi++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          cparam2->cd(ieta+1);
          //g_mPID_m2_pt_ratio[ic][ieta][iphi]->Print();
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->SetStats(0);
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<#eta<%1.1f, %d#pi/6<#phi<%d#pi/6",charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-5)/5.,float(ieta-4)/5.,iphi-6,iphi-5));
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->GetYaxis()->SetTitle("PID Efficiency");
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->SetMarkerStyle(24);
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->Draw("APE");
         
          //m2_PID_pT[icent][ic][ieta][iphi] = new TF1(Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi),m2_PID,0.3,2.5,5);
          nsig_PID_pT[icent][ic][ieta][iphi] = new TF1(Form("nsig_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi),nsig_PID_tanh,0.1,3.5,3);
          nsig_PID_pT[icent][ic][ieta][iphi]->SetParameter(0,0.98);
          nsig_PID_pT[icent][ic][ieta][iphi]->SetParameter(1,20);
          nsig_PID_pT[icent][ic][ieta][iphi]->SetParameter(2,0.0);
          //nsig_PID_pT[icent][ic][ieta][iphi]->SetParameter(3,0.75);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(0,2);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(1,3);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(2,-7);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(3,1.5);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(4,0.5);
          g_mPID_nsig_pt_ratio[ic][ieta][iphi]->Fit(nsig_PID_pT[icent][ic][ieta][iphi],"NMRI");
          nsig_PID_pT[icent][ic][ieta][iphi]->SetLineColor(kRed);
          nsig_PID_pT[icent][ic][ieta][iphi]->Draw("same");
        }
        cparam2->Update();
        cparam2->Print(outputname.c_str());
      }
    }
  }
  cparam2->Print(outputstop.c_str());

  outputname = Form("figures/%s/EfficiencyPID_nsig_1D_eta_ratios.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam2->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ipt = 0; ipt < pt_max; ipt++)
      {
        double pt_low  = 0.1 + 0.2*double(ipt);
        double pt_high = 0.1 + 0.2*double(ipt+1);
        if(ipt == pt_max-1)
        {
          //pt_low = 2.7;
          pt_low = 4.1-(20-pt_max+1)*0.2;
          pt_high = 4.1;
        }
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          cparam2->cd(iphi+1);
          //g_mPID_m2_eta_ratio[ic][ipt][iphi]->Print();
          g_mPID_nsig_eta_ratio[ic][ipt][iphi]->SetStats(0);
          g_mPID_nsig_eta_ratio[ic][ipt][iphi]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<p_{T}<%1.1fGeV/c, %d#pi/6<#phi<%d#pi/6",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,iphi-6,iphi-5));
          g_mPID_nsig_eta_ratio[ic][ipt][iphi]->GetXaxis()->SetTitle("#eta");
          g_mPID_nsig_eta_ratio[ic][ipt][iphi]->GetYaxis()->SetTitle("PID Efficiency");
          g_mPID_nsig_eta_ratio[ic][ipt][iphi]->SetMarkerStyle(24);
          g_mPID_nsig_eta_ratio[ic][ipt][iphi]->Draw("APE");
        }
        cparam2->Update();
        cparam2->Print(outputname.c_str());
      }
    }
  }
  cparam2->Print(outputstop.c_str());

  TCanvas *cparam3 = new TCanvas("c","c",10,10,3000,2500);
  cparam3->Divide(6,5);
  for(int i = 0; i < 30; i++)
  {
    cparam3->cd(i+1);
    cparam3->cd(i+1)->SetLeftMargin(0.15);
    cparam3->cd(i+1)->SetBottomMargin(0.15);
    cparam3->cd(i+1)->SetTicks(1,1);
    cparam3->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_nsig_1D_phi_ratios.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam3->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ipt = 0; ipt < pt_max; ipt++)
      {
        double pt_low  = 0.1 + 0.2*double(ipt);
        double pt_high = 0.1 + 0.2*double(ipt+1);
        if(ipt == pt_max-1)
        {
          pt_low = 4.1-(20-pt_max+1)*0.2;
          pt_high = 4.1;
        }
        for(int ieta = 0; ieta < neta; ieta++)
        {
          cparam3->cd(ieta+1);
          //g_mPID_m2_phi_ratio[ic][ipt][ieta]->Print();
          g_mPID_nsig_phi_ratio[ic][ipt][ieta]->SetStats(0);
          g_mPID_nsig_phi_ratio[ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<p_{T}<%1.1fGeV/c, %1.3f<#eta<%1.3f",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-neta/2)/(neta/2),float(ieta-neta/2)/(neta/2)));
          g_mPID_nsig_phi_ratio[ic][ipt][ieta]->GetXaxis()->SetTitle("#phi");
          g_mPID_nsig_phi_ratio[ic][ipt][ieta]->GetYaxis()->SetTitle("PID Efficiency");
          g_mPID_nsig_phi_ratio[ic][ipt][ieta]->SetMarkerStyle(24);
          g_mPID_nsig_phi_ratio[ic][ipt][ieta]->Draw("APE");
        }
        cparam3->Update();
        cparam3->Print(outputname.c_str());
      }
    }
  }
  cparam3->Print(outputstop.c_str());

  TH1F *parameterhist[10][2][neta][nphi];
  TFile *OutFile = new TFile("nsig_PID_funcparameters_19GeV.root", "RECREATE");
  OutFile->cd();
  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          string histname = Form("nsig_parameters_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi);
          parameterhist[icent][ic][ieta][iphi] = new TH1F(histname.c_str(), histname.c_str(), 3, -0.5, 2.5);
          for(int ipar = 0; ipar < 3; ipar++)
          {
            parameterhist[icent][ic][ieta][iphi]->SetBinContent(ipar+1,nsig_PID_pT[icent][ic][ieta][iphi]->GetParameter(ipar));
          }
          parameterhist[icent][ic][ieta][iphi]->Write();
        }
      }
    }
  }
  OutFile->Close();
 
}
