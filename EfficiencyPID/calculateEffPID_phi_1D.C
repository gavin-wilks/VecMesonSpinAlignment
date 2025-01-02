#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
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

void calculateEffPID_phi_1D(const int energy = 4, const string date = "20240918")
{
  double kaonlimitlow = 0.075;
  double kaonlimitup = 0.1;
  double pionlimitlow = 0.1;
  double pionlimitup = 0.06;

  int pionlowbin = 9;
  int protonlowbin = 11;
 

  //int pt_max = 11;
  //int pt_max_up = 14;
      
  const int pt_lowidx = 1;
  const int pt_highidx = 20; 
  //const int pt_highidx = 2; 
  const int pt_total_bins = 20;
  int pt_rebin_min[pt_total_bins]  = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};//,16};
  int pt_rebin_max[pt_total_bins]  = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  float pt_lower[pt_total_bins] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.5};
  float pt_upper[pt_total_bins]  = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.5,3.5};
  float pt_allbins[pt_total_bins]  = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.5,3.5};

  const int neta = 24;
  float eta_allbins[neta+1];
  for(int ieta = 0; ieta < neta+1; ieta++) eta_allbins[ieta] = float(ieta-neta/2)/float(neta/2);
  for(int ieta = 0; ieta < neta+1; ieta++) cout << eta_allbins[ieta] << endl; // = float(ieta-neta/2)/float(neta/2);

  const int nphi = 24;
  float phi_allbins[nphi+1];
  for(int iphi = 0; iphi < nphi+1; iphi++) phi_allbins[iphi] = TMath::Pi()*float(iphi-nphi/2)/float(nphi/2);
  for(int iphi = 0; iphi < nphi+1; iphi++) cout << phi_allbins[iphi] << endl; // = TMath::Pi()*float(iphi-nphi/2)/float(nphi/2);

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  string InPutFile = Form("data/file_PID_%s_%s_20pt_24eta_24phi.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile = Form("data/Phi_PIDEff_MoreEtaBins.root");
  TFile *File = TFile::Open(InPutFile.c_str());

  TH3F *h_m2_PID[2];
  for(int ic = 0; ic < 2; ic++)
  {
    string histname = Form("h_m2_PID_cent9_charge%d",ic);
    h_m2_PID[ic] = new TH3F(histname.c_str(),histname.c_str(),19,pt_allbins,neta,eta_allbins,nphi,phi_allbins);
  }


  TH1D *h_mPID[2][20][neta][nphi]; // x axis is nsigma_kaon, y axis is m^2
  //TH1D *h_mPID_nsig[2][20][neta][nphi];
  TH1D *h_mPID_m2[2][20][neta][nphi]; 
  TF2 *f_mPID[2][20][neta][nphi];
  //TF1 *f_mPID_nsig[2][20][neta][nphi];
  TF1 *f_mPID_m2[2][20][neta][nphi];
  TFitResultPtr fitResult_mPID_m2[2][20][neta][nphi];
  TF1 *f_1D_mPID[3][2][20][neta][nphi]; // first array runs over particle type: pion, kaon, proton
  //TF1 *f_1D_mPID_nsig[3][2][20][neta][nphi];
  TF1 *f_1D_mPID_m2[3][2][20][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_norm[3][2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_mean[3][2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_sigma[3][2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_pt_ratio[2][neta][nphi];
  TGraphAsymmErrors *g_mPID_m2_eta_ratio[2][20][nphi];
  TGraphAsymmErrors *g_mPID_m2_phi_ratio[2][20][neta];

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
      for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          for(int iphi = 0; iphi < nphi; iphi++)
          {
            cout << "pt = " << ipt << endl;
            for(int ptbin = 0; ptbin < pt_total_bins; ptbin++)
            {
              if(ipt >= pt_rebin_min[ptbin] && ipt < pt_rebin_max[ptbin]) 
              {            
                if(ipt == pt_rebin_min[ptbin])
                {
                  string KEY_InPut = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d_m2",icent,ic,ipt,ieta,iphi);
                  h_mPID_m2[ic][ptbin][ieta][iphi] = (TH1D*)File->Get(KEY_InPut.c_str())->Clone();   
                }
                else
                {
                  string KEY_InPut2 = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d_m2",icent,ic,ipt,ieta,iphi);
                  h_mPID_m2[ic][ptbin][ieta][iphi]->Add((TH1D*)File->Get(KEY_InPut2.c_str())->Clone());   
                }
              } 
            }
          }
        }
      }
    }
  }
  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        cout << "ieta = " << ieta << endl;
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
          {
            ////string KEY_InPut_nsig = Form("h_mPID_nsig_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            //int y_bin_start = h_mPID[ic][ipt][ieta][iphi]->GetXaxis()->FindBin(-2.5);
            //int y_bin_stop  = h_mPID[ic][ipt][ieta][iphi]->GetXaxis()->FindBin(2.5);
            ////h_mPID_nsig[ic][ipt][ieta][iphi] = h_mPID[ic][ipt][ieta][iphi]->ProjectionX(KEY_InPut_nsig.c_str(),y_bin_start,y_bin_stop,"e");
 
            //string KEY_InPut_m2 = Form("h_mPID_m2_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            //h_mPID_m2[ic][ipt][ieta][iphi] = h_mPID[ic][ipt][ieta][iphi]->ProjectionY(KEY_InPut_m2.c_str(),y_bin_start,y_bin_stop,"e");

            string FuncName = Form("f_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
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

            //h_mPID[ic][ipt][ieta]->Fit(f_mPID[ic][ipt][ieta],"NMRI");

            FuncName = Form("f_mPID_m2_cent%d_charge%d_pt%d_eta%d_phi%d",icent,ic,ipt,ieta,iphi);
            if(ipt >= 15) f_mPID_m2[ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),ThreeParticleGaussian,-0.5,1.5,9);
            if(ipt >= protonlowbin && ipt < 15) f_mPID_m2[ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),ThreeParticleGaussian,-0.3,1.2,9);
            if(ipt <  protonlowbin && ipt >= pionlowbin) f_mPID_m2[ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),ThreeParticleGaussian,-0.3,0.5,9);
            if(ipt <  pionlowbin) f_mPID_m2[ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),ThreeParticleGaussian,0.247618-kaonlimitlow,0.247618+kaonlimitup,9);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum());
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(0,3.0,1000000000);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(1,0.01822044);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(1,0.01822044-pionlimitlow,0.01822044+pionlimitup);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(2,0.002+0.006*ipt);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(2,0.001,0.2);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(3,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()/3.0);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(3,3.0,1000000000);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(4,0.247618);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(4,0.247618-kaonlimitlow,0.247618+kaonlimitup);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(5,0.008+0.012*ipt);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(5,0.001,0.2);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(6,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()/2.0);
            if(ipt > 3) f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(6,3.0,1000000000);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(7,0.880354);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(7,0.880354-0.15,0.880354+0.15);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(8,0.044+0.012*ipt);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(8,0.001,0.5);

            //double ptmean = 0.2 + double(ipt)*0.2;
            //if(ipt == pt_max - 1) ptmean = (4.1 - (20 - pt_max_up)*0.2 + (4.1 - (20 - pt_max + 1)*0.2))/2;
            double ptmean = (pt_lower[ipt]+pt_upper[ipt])/2.;//(4.1 - (20 - pt_rebin_max[ipt])*0.2 + (4.1 - (20 - pt_rebin_min[ipt])*0.2))/2;

            cout << "ipt = " << ipt << ", ptmean = " << ptmean << endl; 
            cout << "pt range = [" << pt_lower[ipt] << ", " << pt_upper[ipt] << "}" << endl;

            //if(ipt >= 14)
            //{
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(2,ptmean*0.01,0.2);
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(5,ptmean*0.01,0.2);
            //}

            if(ipt < pionlowbin) 
            {
              f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(0,0.0);
              f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(1,0.0);
              f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(2,0.0);
              //if(ipt < 1) 
              //{
              //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(3,0.0);
              //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(4,0.0);
              //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(5,0.0);
              //}
            }
            if(ipt < protonlowbin) 
            {
              f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(6,0.0);
              f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(7,0.0);
              f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(8,0.0);
              //if(ipt < 1) 
              //{
              //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(3,0.0);
              //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(4,0.0);
              //  f_mPID_m2[ic][ipt][ieta][iphi]->FixParameter(5,0.0);
              //}
            }

            if(ipt >= pt_lowidx+1)  
            {
              //TF1 *tempN[3];           
              //TF1 *tempSigma[3];

              for(int ip = 0; ip < 3; ip++)
              {
                if(ipt < protonlowbin+1 && ip == 2) continue;
                if(ipt < pionlowbin+1  && ip == 0) continue;
                //tempN[ip] = new TF1("tempN","[0]*exp(-[1]*(x-[2]))",0.1,4.1);
                //tempSigma[ip] = new TF1("tempSigma","[0]+[1]*x+[2]*x*x",0.1,4.1);

                double pt, valmean, valsigma, valnorm; 
                g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->GetPoint(ipt-1,pt,valnorm);
                g_mPID_m2_pt_mean[ip][ic][ieta][iphi]->GetPoint(ipt-1,pt,valmean);
                g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->GetPoint(ipt-1,pt,valsigma);

                f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(3*ip  ,valnorm);
                f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(1+3*ip,valmean);
                f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(1+3*ip,valmean-0.03,valmean+0.03);
                f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(2+3*ip,valsigma);
                f_mPID_m2[ic][ipt][ieta][iphi]->SetParLimits(2+3*ip,valsigma-0.01,valsigma+0.075);
             
              }
              //if(ipt >= pionlowbin+1)f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(0,tempN[0]->Eval(ptmean));
              //if(ipt >= pionlowbin+1)f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(2,tempSigma[0]->Eval(ptmean));
              //f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(3,tempN[1]->Eval(ptmean));
              //f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(5,tempSigma[1]->Eval(ptmean));
              //if(ipt >= pionlowbin+1) f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(6,tempN[2]->Eval(ptmean));
              //if(ipt >= pionlowbin+1) f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(8,tempSigma[2]->Eval(ptmean));
            }

            //if(ipt >= 5)  
            //{
            //  TF1 *tempN[3];           
            //  TF1 *tempSigma[3];

            //  for(int ip = 0; ip < 3; ip++)
            //  {
            //    if(ipt < 10 && ip == 2) continue;
            //    if(ipt < 9  && ip == 0) continue;
            //    tempN[ip] = new TF1("tempN","[0]*exp(-[1]*(x-[2]))",0.1,4.1);
            //    tempSigma[ip] = new TF1("tempSigma","[0]+[1]*x+[2]*x*x",0.1,4.1);

            //    double pt, val; 
            //    if(ip == 0) g_mPID_m2_pt_mean[0][ic][ieta][iphi]->GetPoint(ipt,pt,val);
            //    if(ip == 1) g_mPID_m2_pt_mean[1][ic][ieta][iphi]->GetPoint(ipt,pt,val);
            //    if(ip == 2) g_mPID_m2_pt_mean[2][ic][ieta][iphi]->GetPoint(ipt,pt,val);
            //    //tempN[ip]->SetParameter(0,val);

            //    g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->Fit(tempN[ip],"QNRI");
            //    g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->Fit(tempSigma[ip],"QNRI");

            //    f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(1+3*ip,val);
            // 
            //  }
            //  if(ipt >= 9)f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(0,tempN[0]->Eval(ptmean));
            //  if(ipt >= 9)f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(2,tempSigma[0]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(3,tempN[1]->Eval(ptmean));
            //  f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(5,tempSigma[1]->Eval(ptmean));
            //  if(ipt >= 10) f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(6,tempN[2]->Eval(ptmean));
            //  if(ipt >= 10) f_mPID_m2[ic][ipt][ieta][iphi]->SetParameter(8,tempSigma[2]->Eval(ptmean));
            //}
        

            fitResult_mPID_m2[ic][ipt][ieta][iphi] = h_mPID_m2[ic][ipt][ieta][iphi]->Fit(f_mPID_m2[ic][ipt][ieta][iphi],"QNRIBS");


            for(int ip = 0; ip < 3; ip++)
            {

              g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->SetPoint(ipt,ptmean,  f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(3*ip));
              g_mPID_m2_pt_norm[ip][ic][ieta][iphi]->SetPointError(ipt,0.0,0.0, f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(3*ip),  f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(3*ip));
              g_mPID_m2_pt_mean[ip][ic][ieta][iphi]->SetPoint(ipt,ptmean,  f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(1+3*ip));
              g_mPID_m2_pt_mean[ip][ic][ieta][iphi]->SetPointError(ipt,0.0,0.0, f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(1+3*ip),f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(1+3*ip));
              g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->SetPoint(ipt,ptmean, f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(2+3*ip));
              g_mPID_m2_pt_sigma[ip][ic][ieta][iphi]->SetPointError(ipt,0.0,0.0,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(2+3*ip),f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(2+3*ip));
            }
            //FuncName = Form("f_mPID_nsig_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
            //f_mPID_nsig[ic][ipt][ieta] = new TF1(FuncName.c_str(),ThreeParticleGaussian,-10.0,10.0,9);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(0,h_mPID_nsig[ic][ipt][ieta]->GetMaximum());
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(1,0.01822044);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParLimits(1,0.01822044-0.05,0.01822044+0.05);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(2,0.1);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(3,h_mPID_nsig[ic][ipt][ieta]->GetMaximum()/3.0);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(4,0.247618);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParLimits(4,0.247618-0.05,0.247618+0.05);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(5,0.1);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(6,h_mPID_nsig[ic][ipt][ieta]->GetMaximum()/2.0);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(7,0.880354);
            ////f_mPID_nsig[ic][ipt][ieta]->SetParLimits(7,0.880354-0.05,0.880354+0.05);
            //f_mPID_nsig[ic][ipt][ieta]->SetParameter(8,0.1);

            //h_mPID_nsig[ic][ipt][ieta]->Fit(f_mPID_nsig[ic][ipt][ieta],"NMRI");
          }
        }
      }
      
    }
  }

  string particlename[3] = {"pion","kaon","proton"};
  string charge[2] = {"K-","K+"};
  int cent_low[10]  = {70,60,50,40,30,20,10,5,0,20};
  int cent_high[10] = {80,70,60,50,40,30,20,10,5,60};

  TCanvas *cparam = new TCanvas("c","c",10,10,2400,1600);
  cparam->Divide(6,4);
  for(int i = 0; i < 24; i++)
  {
    cparam->cd(i+1);
    cparam->cd(i+1)->SetLeftMargin(0.15);
    cparam->cd(i+1)->SetBottomMargin(0.15);
    cparam->cd(i+1)->SetTicks(1,1);
    cparam->cd(i+1)->SetGrid(0,0);
  }

  string outputname = Form("figures/%s/EfficiencyPID_m2_1D_means_1DHistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
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
  //          g_mPID_m2_pt_mean[ipart][ic][ieta][iphi]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-neta/2)/float(neta/2),float(ieta-neta/2)/float(neta/2)));
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

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_sigmas_1Dhistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
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
  //          g_mPID_m2_pt_sigma[ipart][ic][ieta][iphi]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-neta/2)/float(neta/2),float(ieta-neta/2)/float(neta/2)));
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

  //outputname = Form("figures/%s/EfficiencyPID_m2_1D_norms_1Dhistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
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
  //          g_mPID_m2_pt_norm[ipart][ic][ieta][iphi]->SetTitle(Form("%s, %s, Cent=%d-%d, %1.1f<#eta<%1.1f",particlename[ipart].c_str(),charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-neta/2)/float(neta/2),float(ieta-neta/2)/float(neta/2)));
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


  TCanvas *c = new TCanvas("c","c",10,10,2400,1600);
  c->Divide(6,4);
  for(int i = 0; i < 12; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_nsigma_1D_1DHistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  //c->Print(outputstart.c_str());

  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = pt_lowidx; ipt < pt_max; ipt++)
  //    {
  //      for(int ieta = 0; ieta < 20; ieta++)
  //      {
  //        c->cd(ieta+1);
  //        h_mPID_nsig[ic][ipt][ieta]->SetStats(0);
  //        double pt_low  = 0.1 + 0.2*double(ipt);
  //        double pt_high = 0.1 + 0.2*double(ipt+1);
  //        h_mPID_nsig[ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high));
  //        h_mPID_nsig[ic][ipt][ieta]->GetXaxis()->SetTitle("n#sigma_{K}");
  //        h_mPID_nsig[ic][ipt][ieta]->GetYaxis()->SetTitle("Number of Tracks");
  //        h_mPID_nsig[ic][ipt][ieta]->Draw("pE");
  //      }
  //      c->Update();
  //      c->Print(outputname.c_str());
  //    }
  //  }
  //}
  //c->Print(outputstop.c_str());


  outputname = Form("figures/%s/EfficiencyPID_m2_1D_1DHistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  int particlecolor[3] = {kOrange+7,kBlue,kGray+2};

  c->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          g_mPID_m2_pt_ratio[ic][ieta][iphi] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
      {
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          g_mPID_m2_eta_ratio[ic][ipt][iphi] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          g_mPID_m2_phi_ratio[ic][ipt][ieta] = new TGraphAsymmErrors();
        }
      }
      for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
      {
        for(int ieta = 0; ieta < neta; ieta++)
        {
          for(int iphi = 0; iphi < nphi; iphi++)
          {
            c->cd(iphi+1);
            //c->cd(iphi+1)->SetLogy();
            h_mPID_m2[ic][ipt][ieta][iphi]->SetStats(0);
            //double pt_low  = 0.1 + 0.2*double(pt_rebin_min[ipt]);
            //double pt_high = 0.1 + 0.2*double(pt_rebin_max[ipt]);
            double pt_low  = pt_lower[ipt];//0.1 + 0.2*double(pt_rebin_min[ipt]);
            double pt_high = pt_upper[ipt];//0.1 + 0.2*double(pt_rebin_max[ipt]);
            //double pt_low  = 0.1 + 0.2*double(ipt);
            //double pt_high = 0.1 + 0.2*double(ipt+1);
            //if(ipt == pt_max - 1)
            //{
            //  ///pt_low = 2.7;
            //  pt_low = 4.1 - (20 - pt_max + 1) * 0.2;
            //  pt_high = 4.1 - (20 - pt_max_up) * 0.2;
            //}
            double etamean = (double(ieta)-double(neta/2)+0.5)/double(neta/2);
            double phimean = TMath::Pi()*(double(iphi)-double(nphi/2)+0.5)/double(nphi/2);

            h_mPID_m2[ic][ipt][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c, %1.1f<#eta<%1.1f, %d#pi/%d<#phi<%d#pi/%d",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-neta/2)/float(neta/2),float(ieta-neta/2)/float(neta/2),iphi-nphi/2,nphi/2,iphi-nphi/2+1,nphi/2));
            h_mPID_m2[ic][ipt][ieta][iphi]->GetXaxis()->SetTitle("m^{2} (GeV/c^{2})");
            h_mPID_m2[ic][ipt][ieta][iphi]->GetYaxis()->SetTitle("Number of Tracks");
            h_mPID_m2[ic][ipt][ieta][iphi]->Draw("EX0");

            f_mPID_m2[ic][ipt][ieta][iphi]->SetLineColor(kMagenta);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetNpx(1000);
            f_mPID_m2[ic][ipt][ieta][iphi]->SetLineWidth(1);
            f_mPID_m2[ic][ipt][ieta][iphi]->Draw("same");
      
            for(int ipart = 0; ipart < 3; ipart++)
            {
              string FuncName = Form("f_mPID_m2_%s_cent%d_charge%d_pt%d_eta%d_phi%d",particlename[ipart].c_str(),icent,ic,ipt,ieta,iphi);
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi] = new TF1(FuncName.c_str(),SingleParticleGaussian,-0.5,1.5,3);
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetNpx(1000);
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParameter(0,f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(0+ipart*3));
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParameter(1,f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(1+ipart*3));
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParameter(2,f_mPID_m2[ic][ipt][ieta][iphi]->GetParameter(2+ipart*3));
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParError(0,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(0+ipart*3));
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParError(1,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(1+ipart*3));
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetParError(2,f_mPID_m2[ic][ipt][ieta][iphi]->GetParError(2+ipart*3));

              //if(ipart == 1 && ipt > 0)
              if(ipart == 1)
              {
                double params[3] = {0.0};
                TMatrixDSym covMat(3);
                params[0] = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetParams()[3];
                params[1] = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetParams()[4];
                params[2] = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetParams()[5];
                covMat(0,0) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(3,3);
                covMat(0,1) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(3,4);
                covMat(0,2) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(3,5);
                covMat(1,0) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(4,3);
                covMat(1,1) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(4,4);
                covMat(1,2) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(4,5);
                covMat(2,0) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(5,3);
                covMat(2,1) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(5,4);
                covMat(2,2) = fitResult_mPID_m2[ic][ipt][ieta][iphi]->GetCovarianceMatrix()(5,5);

                float bin_width = h_mPID_m2[ic][ipt][ieta][iphi]->GetBinWidth(1);

                double kaons_selected  = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->Integral(0.16,0.36)/bin_width;   
                double dkaons_selected = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->IntegralError(0.16,0.36,params,covMat.GetMatrixArray())/bin_width;   
                double kaons_total     = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->Integral(-0.5,0.75)/bin_width;   
                double dkaons_total    = f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->IntegralError(-0.5,0.75,params,covMat.GetMatrixArray())/bin_width;
         
                //if(kaons_total == 0.0) continue; 

                double covariance = dkaons_selected * dkaons_selected;

                //cout << "kaons_selected = " << kaons_selected << " +/- " << dkaons_selected << endl;
                //cout << "kaons_total    = " << kaons_total << " +/- " << dkaons_total << endl;
                
                double R = kaons_selected/kaons_total;
                if(R > 1 + 1e-4 || R < 0) 
                {
                  cout << "kaons_selected = " << kaons_selected << " +/- " << dkaons_selected << endl;
                  cout << "kaons_total    = " << kaons_total << " +/- " << dkaons_total << endl;
                  cout << "R = " << R << endl;
                  continue;               
                }
      
                double relError_selected = dkaons_selected/kaons_selected;   
                double relError_total    = dkaons_total   /kaons_total;   
                double correlatedTerm = -2.0*covariance/(kaons_selected*kaons_total);

                //cout << "relError_selected = " << relError_selected << endl;
                //cout << "relError_total    = " << relError_total << endl;
                //cout << "correlatedTerm    = " << correlatedTerm << endl;
                double sqrtinnerds = relError_selected * relError_selected + relError_total * relError_total + correlatedTerm;
                //cout << "sqrtinnerds = " << sqrtinnerds << endl;
                double deltaR = R * sqrt(sqrtinnerds);
                if(sqrtinnerds < 1e-10 && sqrtinnerds > -1e-10) deltaR = 0.0;
                //if(deltaR > 1) continue;
                if(deltaR > 1 || std::isnan(deltaR)) {

                  cout << "relError_selected = " << relError_selected << endl;
                  cout << "relError_total    = " << relError_total << endl;
                  cout << "correlatedTerm    = " << correlatedTerm << endl;
                  cout << "sqrtinnerds = " << sqrtinnerds << endl;
                
                  deltaR = 0;      
                  cout << "Filling plots with " << R << " +/- " << deltaR << endl;
        
                }
            
                h_m2_PID[ic]->SetBinContent(ipt-pt_lowidx+1,ieta+1,iphi+1,R);   
                h_m2_PID[ic]->SetBinError(ipt-pt_lowidx+1,ieta+1,iphi+1,deltaR);   

                g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetPoint(ipt,(pt_low+pt_high)/2.0,R);
                g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetPointError(ipt,0.0,0.0,deltaR,deltaR);
                g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetPoint(ieta,etamean,R);
                g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetPointError(ieta,0.0,0.0,deltaR,deltaR);
                g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetPoint(iphi,phimean,R);
                g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetPointError(iphi,0.0,0.0,deltaR,deltaR);
              }

	      f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetLineColor(particlecolor[ipart]);
	      f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetLineWidth(1);
	      f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->SetLineStyle(1);
              f_1D_mPID_m2[ipart][ic][ipt][ieta][iphi]->Draw("same");

              PlotLine(0.01822044-pionlimitlow,0.01822044-pionlimitlow,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kBlack,1,2);
              PlotLine(0.01822044+pionlimitup,0.01822044+pionlimitup,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kBlack,1,2);
              PlotLine(0.247618-kaonlimitlow,0.247618-kaonlimitlow,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
              PlotLine(0.247618+kaonlimitup,0.247618+kaonlimitup,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
              PlotLine(0.880354-0.15,0.880354-0.15,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
              PlotLine(0.880354+0.15,0.880354+0.15,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum()    ,kBlack,1,2);
 
              PlotLine(0.16,0.16,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kGreen,1,2);
              PlotLine(0.36,0.36,0,h_mPID_m2[ic][ipt][ieta][iphi]->GetMaximum(),kGreen,1,2);
            }
          } 
          c->Update();
          c->Print(outputname.c_str());
        }
      }
    }
  }
  c->Print(outputstop.c_str());

  cout << "Finished with plots before ratios" << endl;

  TCanvas *cparam2 = new TCanvas("c","c",10,10,2400,1600);
  cparam2->Divide(6,4);
  for(int i = 0; i < 12; i++)
  {
    cparam2->cd(i+1);
    cparam2->cd(i+1)->SetLeftMargin(0.15);
    cparam2->cd(i+1)->SetBottomMargin(0.15);
    cparam2->cd(i+1)->SetTicks(1,1);
    cparam2->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_m2_1D_pt_ratios_1Dhistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam2->Print(outputstart.c_str());

  TF1 *m2_PID_pT[10][2][neta][nphi];

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          cparam2->cd(iphi+1);
          //g_mPID_m2_pt_ratio[ic][ieta][iphi]->Print();
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetStats(0);
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<#eta<%1.1f, %d#pi/%d<#phi<%d#pi/%d",charge[ic].c_str(),cent_low[icent],cent_high[icent],float(ieta-neta/2)/float(neta/2),float(ieta-neta/2)/float(neta/2),iphi-nphi/2,nphi/2,iphi-nphi/2+1,nphi/2));
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->GetYaxis()->SetTitle("PID Efficiency");
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->SetMarkerStyle(24);
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->Draw("APE");
         
          //m2_PID_pT[icent][ic][ieta][iphi] = new TF1(Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi),m2_PID,0.3,2.5,5);
          m2_PID_pT[icent][ic][ieta][iphi] = new TF1(Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi),m2_PID_tanh,1.1,2.7,4);
          m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(0,0.25);
          m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(1,-2);
          m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(2,1.75);
          m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(3,0.75);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(0,2);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(1,3);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(2,-7);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(3,1.5);
          //m2_PID_pT[icent][ic][ieta][iphi]->SetParameter(4,0.5);
          g_mPID_m2_pt_ratio[ic][ieta][iphi]->Fit(m2_PID_pT[icent][ic][ieta][iphi],"NRI");
          m2_PID_pT[icent][ic][ieta][iphi]->SetLineColor(kRed);
          m2_PID_pT[icent][ic][ieta][iphi]->Draw("same");
        }
        cparam2->Update();
        cparam2->Print(outputname.c_str());
      }
    }
  }
  cparam2->Print(outputstop.c_str());

  outputname = Form("figures/%s/EfficiencyPID_m2_1D_eta_ratios_1Dhistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam2->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
      {
        double pt_low  = pt_lower[ipt];//0.1 + 0.2*double(pt_rebin_min[ipt]);
        double pt_high = pt_upper[ipt];//0.1 + 0.2*double(pt_rebin_max[ipt]);
        //double pt_low  = 0.1 + 0.2*double(pt_rebin_min[ipt]);
        //double pt_high = 0.1 + 0.2*double(pt_rebin_max[ipt]);
        //double pt_low  = 0.1 + 0.2*double(ipt);
        //double pt_high = 0.1 + 0.2*double(ipt+1);
        //if(ipt == pt_max - 1)
        //{
        //  pt_low = 4.1 - (20 - pt_max +1) * 0.2;
        //  //pt_low = 2.7;
        //  pt_high = 4.1 - (20 - pt_max_up) * 0.2;
        //}
        for(int iphi = 0; iphi < nphi; iphi++)
        {
          cparam2->cd(iphi+1);
          //g_mPID_m2_eta_ratio[ic][ipt][iphi]->Print();
          g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetStats(0);
          g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<p_{T}<%1.1fGeV/c, %d#pi/%d<#phi<%d#pi/%d",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,iphi-nphi/2,nphi/2,iphi-nphi/2+1,nphi/2));
          g_mPID_m2_eta_ratio[ic][ipt][iphi]->GetXaxis()->SetTitle("#eta");
          g_mPID_m2_eta_ratio[ic][ipt][iphi]->GetYaxis()->SetTitle("PID Efficiency");
          g_mPID_m2_eta_ratio[ic][ipt][iphi]->SetMarkerStyle(24);
          g_mPID_m2_eta_ratio[ic][ipt][iphi]->Draw("APE");
        }
        cparam2->Update();
        cparam2->Print(outputname.c_str());
      }
    }
  }
  cparam2->Print(outputstop.c_str());

  TCanvas *cparam3 = new TCanvas("c","c",10,10,2400,1600);
  cparam3->Divide(6,4);
  for(int i = 0; i < neta; i++)
  {
    cparam3->cd(i+1);
    cparam3->cd(i+1)->SetLeftMargin(0.15);
    cparam3->cd(i+1)->SetBottomMargin(0.15);
    cparam3->cd(i+1)->SetTicks(1,1);
    cparam3->cd(i+1)->SetGrid(0,0);
  }

  outputname = Form("figures/%s/EfficiencyPID_m2_1D_phi_ratios_1DHistograms.pdf",vmsa::mBeamEnergy[energy].c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop  = Form("%s]",outputname.c_str());

  cparam3->Print(outputstart.c_str());

  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ipt = pt_lowidx; ipt < pt_highidx; ipt++)
      {
        double pt_low  = pt_lower[ipt];//0.1 + 0.2*double(pt_rebin_min[ipt]);
        double pt_high = pt_upper[ipt];//0.1 + 0.2*double(pt_rebin_max[ipt]);
        //double pt_low  = 0.1 + 0.2*double(pt_rebin_min[ipt]);
        //double pt_high = 0.1 + 0.2*double(pt_rebin_max[ipt]);
        //double pt_low  = 0.1 + 0.2*double(ipt);
        //double pt_high = 0.1 + 0.2*double(ipt+1);
        //if(ipt == pt_max - 1 )
        //{
        //  //pt_low = 2.7;
        //  pt_low = 4.1 - (20 - pt_max + 1) * 0.2;
        //  pt_high = 4.1 - (20 - pt_max_up) * 0.2;
        //}
        for(int ieta = 0; ieta < neta; ieta++)
        {
          cparam3->cd(ieta+1);
          //g_mPID_m2_phi_ratio[ic][ipt][ieta]->Print();
          g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetStats(0);
          g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, %1.1f<p_{T}<%1.1fGeV/c, %1.1f<#eta<%1.1f",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high,float(ieta-neta/2)/float(neta/2),float(ieta-neta/2)/float(neta/2)));
          g_mPID_m2_phi_ratio[ic][ipt][ieta]->GetXaxis()->SetTitle("#phi");
          g_mPID_m2_phi_ratio[ic][ipt][ieta]->GetYaxis()->SetTitle("PID Efficiency");
          g_mPID_m2_phi_ratio[ic][ipt][ieta]->SetMarkerStyle(24);
          g_mPID_m2_phi_ratio[ic][ipt][ieta]->Draw("APE");
        }
        cparam3->Update();
        cparam3->Print(outputname.c_str());
      }
    }
  }
  cparam3->Print(outputstop.c_str());

  //outputname = Form("figures/%s/EfficiencyPID_2D.pdf",vmsa::mBeamEnergy[energy].c_str());
  //outputstart = Form("%s[",outputname.c_str());
  //outputstop  = Form("%s]",outputname.c_str());

  //c->Print(outputstart.c_str());

  //for(int icent = 0; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ipt = pt_lowidx; ipt < pt_max; ipt++)
  //    {
  //      for(int ieta = 0; ieta < 20; ieta++)
  //      {
  //        c->cd(ieta+1);
  //        c->cd(ieta+1)->SetRightMargin(0.15);
  //        h_mPID[ic][ipt][ieta]->SetStats(0);
  //        double pt_low  = 0.1 + 0.2*double(ipt);
  //        double pt_high = 0.1 + 0.2*double(ipt+1);
  //        h_mPID[ic][ipt][ieta]->SetTitle(Form("%s, Cent=%d-%d, p_{T}=[%1.1f,%1.1f} GeV/c",charge[ic].c_str(),cent_low[icent],cent_high[icent],pt_low,pt_high));
  //        h_mPID[ic][ipt][ieta]->GetXaxis()->SetTitle("n#sigma_{K}");
  //        h_mPID[ic][ipt][ieta]->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})");
  //        h_mPID[ic][ipt][ieta]->GetZaxis()->SetTitle("Number of Tracks");
  //        h_mPID[ic][ipt][ieta]->Draw("colz");
  //      }
  //      c->Update();
  //      c->Print(outputname.c_str());
  //    }
  //  }
  //}
  //c->Print(outputstop.c_str());

  //int npar = 4;
  //TH1F *parameterhist[10][2][neta][nphi];
  //TFile *OutFile = new TFile("m2_PID_funcparameters_19GeV.root", "RECREATE");
  //OutFile->cd();
  //for(int icent = 9; icent < 10; icent++)
  //{
  //  for(int ic = 0; ic < 2; ic++)
  //  {
  //    for(int ieta = 0; ieta < neta; ieta++)
  //    {
  //      for(int iphi = 0; iphi < nphi; iphi++)
  //      {
  //        string histname = Form("m2_parameters_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi);
  //        parameterhist[icent][ic][ieta][iphi] = new TH1F(histname.c_str(), histname.c_str(), npar, -0.5, double(npar)-0.5);
  //        for(int ipar = 0; ipar < npar; ipar++)
  //        {
  //          parameterhist[icent][ic][ieta][iphi]->SetBinContent(ipar+1,m2_PID_pT[icent][ic][ieta][iphi]->GetParameter(ipar));
  //        }
  //        parameterhist[icent][ic][ieta][iphi]->Write();
  //      }
  //    }
  //  }
  //}
  //OutFile->Close();
  TFile *OutFile = new TFile("m2_PID_hist_19GeV.root", "RECREATE");
  OutFile->cd();
  for(int icent = 9; icent < 10; icent++)
  {
    for(int ic = 0; ic < 2; ic++)
    {
      h_m2_PID[ic]->Write();
    }
  }
  OutFile->Close();
}
