#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
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
//#ifdef MAKECINT
//#pragma link C++ class std::map<std::string,TH1F*>+;
//#endif

#ifndef _PlotQA_
#define _PlotQA_  1
#endif
#ifndef _SaveQA_
#define _SaveQA_  0
#endif

using namespace std;

void Rebin(TH1F* hRebinned)
{
  for (int i = 1; i <= hRebinned->GetNbinsX(); ++i) 
  {
    double binWidth = hRebinned->GetBinWidth(i);
    double binContent = hRebinned->GetBinContent(i);
    double binError = hRebinned->GetBinError(i);

    hRebinned->SetBinContent(i, binContent / binWidth);
    hRebinned->SetBinError(i, binError / binWidth);
  }
}

void compare_CWR_GAVIN(int energy = 5, int pid = 0, int year = 0, string date = "20240424", bool random3D = false, int order = 2, string etamode = "eta1_eta1")
{
  TGaxis::SetMaxDigits(4);

  string inputfileCWR   = Form("../output/AuAu%s/%s/output_2718_pTheRct.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfileGAVIN = Form("../output/AuAu%s/%s/AllRawHists_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode.c_str());
  TFile *File_CWR = TFile::Open(inputfileCWR.c_str());
  TFile *File_GAVIN = TFile::Open(inputfileGAVIN.c_str());

  TH1FMap h_mGAVIN_SE, h_mCWR_SE;
  TH1FMap h_mGAVIN_ME, h_mCWR_ME;

  int const pt_rebin_local = 6; // maximum pt binning
  float const pt_low_local[pt_rebin_local] = {0.0,1.2,1.8,2.4,3.0,4.2};
  float const pt_up_local[pt_rebin_local]  = {1.2,1.8,2.4,3.0,4.2,5.4};
  int const pt_rebin_start_local[pt_rebin_local] = {0,5,8,11,14,17};
  int const pt_rebin_stop_local[pt_rebin_local]  = {4,7,10,13,16,19};

  // Load Gavin's Plots and Rebin pT 
  for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
  {
    string KEY_SE_Int = Form("Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",i_theta);
    string KEY_ME_Int = Form("Centrality_9_CosThetaStar_%d_ME_Rebinned_GAVIN",i_theta);
    for(int pt_bin = 0; pt_bin < pt_rebin_local; pt_bin++) // pt loop
    {
      string KEY_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",pt_bin,i_theta);
      string KEY_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_ME_Rebinned_GAVIN",pt_bin,i_theta);
      for(int i_pt = pt_rebin_start_local[pt_bin]; i_pt <= pt_rebin_stop_local[pt_bin]; i_pt++)
      {
        string KEY_GAVIN_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE",i_pt,i_theta);
        string KEY_GAVIN_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_ME",i_pt,i_theta);

        h_mGAVIN_SE[KEY_GAVIN_SE] = (TH1F*) File_GAVIN->Get(KEY_GAVIN_SE.c_str());
        h_mGAVIN_ME[KEY_GAVIN_ME] = (TH1F*) File_GAVIN->Get(KEY_GAVIN_ME.c_str());
        h_mGAVIN_SE[KEY_GAVIN_SE]->Rebin(2);
        h_mGAVIN_ME[KEY_GAVIN_ME]->Rebin(2);
        h_mGAVIN_SE[KEY_GAVIN_SE]->Scale(1.89/2.443373);
        h_mGAVIN_ME[KEY_GAVIN_ME]->Scale(1.89/2.443373);
        
        //Rebin(h_mGAVIN_SE[KEY_GAVIN_SE]);
        //Rebin(h_mGAVIN_ME[KEY_GAVIN_ME]);

        //cout << "nbins GAVIN = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetNbinsX() << ", bin width = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinWidth(1) << ", low edge = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinLowEdge(1) << ", high edge = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinLowEdge(h_mGAVIN_SE[KEY_GAVIN_SE]->GetNbinsX()) + h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinWidth(1) << endl;        

        if(i_pt == pt_rebin_start_local[pt_bin])
        {
          h_mGAVIN_SE[KEY_SE] = (TH1F*) h_mGAVIN_SE[KEY_GAVIN_SE]->Clone(KEY_SE.c_str()); 
          h_mGAVIN_ME[KEY_ME] = (TH1F*) h_mGAVIN_ME[KEY_GAVIN_ME]->Clone(KEY_ME.c_str()); 
        }
        else
        {
          h_mGAVIN_SE[KEY_SE]->Add(h_mGAVIN_SE[KEY_GAVIN_SE],1.0); 
          h_mGAVIN_ME[KEY_ME]->Add(h_mGAVIN_ME[KEY_GAVIN_ME],1.0); 
        }
        // pt 1.2-5.4
        if(i_pt == pt_rebin_start_local[1])
        {
          h_mGAVIN_SE[KEY_SE_Int] = (TH1F*) h_mGAVIN_SE[KEY_GAVIN_SE]->Clone(KEY_SE_Int.c_str()); 
          h_mGAVIN_ME[KEY_ME_Int] = (TH1F*) h_mGAVIN_ME[KEY_GAVIN_ME]->Clone(KEY_ME_Int.c_str()); 
        }
        else if(i_pt > pt_rebin_start_local[1] && i_pt <= pt_rebin_stop_local[5])
        {
          h_mGAVIN_SE[KEY_SE_Int]->Add(h_mGAVIN_SE[KEY_GAVIN_SE],1.0); 
          h_mGAVIN_ME[KEY_ME_Int]->Add(h_mGAVIN_ME[KEY_GAVIN_ME],1.0); 
        }
      }
    }
  }

  TH1FMap h_mRatio_SE, h_mRatio_ME;

  int const pt_rebin_start_local_CWR[pt_rebin_local] = {0,1,2,3,4,6};
  int const pt_rebin_stop_local_CWR[pt_rebin_local]  = {0,1,2,3,5,7};
  

  // Load CWR's Plots and Rebin pT 
  for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
  {
    string KEY_SE_Int = Form("Centrality_9_CosThetaStar_%d_SE_Rebinned_CWR",i_theta);
    string KEY_ME_Int = Form("Centrality_9_CosThetaStar_%d_ME_Rebinned_CWR",i_theta);
    string KEY_GAVIN_SE_Int = Form("Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",i_theta);
    string KEY_GAVIN_ME_Int = Form("Centrality_9_CosThetaStar_%d_ME_Rebinned_GAVIN",i_theta);
    cout << "KEY_SE_Int = " << KEY_SE_Int << endl;
    for(int pt_bin = 0; pt_bin < pt_rebin_local; pt_bin++) // pt loop
    {
      string KEY_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_CWR",pt_bin,i_theta);
      string KEY_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_ME_Rebinned_CWR",pt_bin,i_theta);
      string KEY_GAVIN_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",pt_bin,i_theta);
      string KEY_GAVIN_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_ME_Rebinned_GAVIN",pt_bin,i_theta);
      for(int i_pt = pt_rebin_start_local_CWR[pt_bin]; i_pt <= pt_rebin_stop_local_CWR[pt_bin]; i_pt++)
      {
        for(int i_cent = 3; i_cent <= 6; i_cent++)
        {
          string KEY_CWR_SE = Form("DefRelPairMassOSCent%dPt%dCosT%d",i_cent,i_pt,i_theta);
          string KEY_CWR_ME = Form("DefMixPairMassOSCent%dPt%dCosT%d",i_cent,i_pt,i_theta);

          TH1F *SE_temp = (TH1F*) File_CWR->Get(KEY_CWR_SE.c_str());
          TH1F *ME_temp = (TH1F*) File_CWR->Get(KEY_CWR_ME.c_str());

          string KEY_temp_SE = Form("pt_%d_Centrality_%d_CosThetaStar_%d_SE_CWR",i_pt,i_cent,i_theta);
          string KEY_temp_ME = Form("pt_%d_Centrality_%d_CosThetaStar_%d_ME_CWR",i_pt,i_cent,i_theta);

          h_mCWR_SE[KEY_CWR_SE] = new TH1F(KEY_temp_SE.c_str(),KEY_temp_SE.c_str(),350,0.98,1.05);
          h_mCWR_ME[KEY_CWR_ME] = new TH1F(KEY_temp_ME.c_str(),KEY_temp_ME.c_str(),350,0.98,1.05);

          for(int ibin = 401; ibin <= 750; ibin++)
          {
            double valSE = SE_temp->GetBinContent(ibin);
            double errSE = SE_temp->GetBinError(ibin);
            //cout << "low edge = " << SE_temp->GetBinLowEdge(ibin) << endl;
            //cout << "valSE = " << valSE << " +/- " << errSE << endl;
            double valME = ME_temp->GetBinContent(ibin);
            double errME = ME_temp->GetBinError(ibin);

            h_mCWR_SE[KEY_CWR_SE]->SetBinContent(ibin-400,valSE);
            h_mCWR_SE[KEY_CWR_SE]->SetBinError(ibin-400,errSE);
            h_mCWR_ME[KEY_CWR_ME]->SetBinContent(ibin-400,valME);
            h_mCWR_ME[KEY_CWR_ME]->SetBinError(ibin-400,errME);
            double valSEnew = h_mCWR_SE[KEY_CWR_SE]->GetBinContent(ibin-400);
            double errSEnew = h_mCWR_SE[KEY_CWR_SE]->GetBinError(ibin-400);
            //cout << "valSEnew = " << valSEnew << " +/- " << errSEnew << endl;
          }  
          h_mCWR_SE[KEY_CWR_SE]->Rebin(7);
          h_mCWR_ME[KEY_CWR_ME]->Rebin(7);
          //Rebin(h_mCWR_SE[KEY_CWR_SE]);
          //Rebin(h_mCWR_ME[KEY_CWR_ME]);
            
          //cout << "nbins CWR = " << h_mCWR_SE[KEY_CWR_SE]->GetNbinsX() << ", bin width = " << h_mCWR_SE[KEY_CWR_SE]->GetBinWidth(1) << ", low edge = " << h_mCWR_SE[KEY_CWR_SE]->GetBinLowEdge(1) << ", high edge = " << h_mCWR_SE[KEY_CWR_SE]->GetBinLowEdge(h_mCWR_SE[KEY_CWR_SE]->GetNbinsX()) + h_mCWR_SE[KEY_CWR_SE]->GetBinWidth(1) << endl;        

          if(i_pt == pt_rebin_start_local_CWR[pt_bin] && i_cent == 3)
          {
            h_mCWR_SE[KEY_SE] = (TH1F*) h_mCWR_SE[KEY_CWR_SE]->Clone(KEY_SE.c_str()); 
            h_mCWR_ME[KEY_ME] = (TH1F*) h_mCWR_ME[KEY_CWR_ME]->Clone(KEY_ME.c_str()); 
          }
          else
          {
            h_mCWR_SE[KEY_SE]->Add(h_mCWR_SE[KEY_CWR_SE],1.0); 
            h_mCWR_ME[KEY_ME]->Add(h_mCWR_ME[KEY_CWR_ME],1.0); 
          }
          // pt 1.2-5.4
          if(i_pt == pt_rebin_start_local_CWR[1] && i_cent == 3)
          {
            cout << "Creating plot for ipt = " << i_pt << ", icent = " << i_cent << endl;
            cout << "# of entries = " << h_mCWR_SE[KEY_CWR_SE]->Integral(0,-1) << endl;
            h_mCWR_SE[KEY_SE_Int] = (TH1F*) h_mCWR_SE[KEY_CWR_SE]->Clone(KEY_SE_Int.c_str()); 
            h_mCWR_ME[KEY_ME_Int] = (TH1F*) h_mCWR_ME[KEY_CWR_ME]->Clone(KEY_ME_Int.c_str()); 
          }
          else if((i_pt == pt_rebin_start_local_CWR[1] && i_cent > 3) || (i_pt > pt_rebin_start_local_CWR[1] && i_cent >= 3))
          {
            cout << "Adding plot for ipt = " << i_pt << ", icent = " << i_cent << endl;
            cout << "# of entries = " << h_mCWR_SE[KEY_CWR_SE]->Integral(0,-1) << endl;
            h_mCWR_SE[KEY_SE_Int]->Add(h_mCWR_SE[KEY_CWR_SE],1.0); 
            h_mCWR_ME[KEY_ME_Int]->Add(h_mCWR_ME[KEY_CWR_ME],1.0); 
          }
        }
      }
      string KEY_ratio_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Ratio",pt_bin,i_theta);
      string KEY_ratio_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_ME_Ratio",pt_bin,i_theta);
      h_mRatio_SE[KEY_ratio_SE] = (TH1F*) h_mCWR_SE[KEY_SE]->Clone(KEY_ratio_SE.c_str());
      h_mRatio_SE[KEY_ratio_SE]->Divide(h_mGAVIN_SE[KEY_GAVIN_SE]);
      //h_mRatio_SE[KEY_ratio_SE]->Divide(h_mCWR_SE[KEY_SE],h_mGAVIN_SE[KEY_GAVIN_SE],1,1,"B");
      //cout << "nbins CWR = " << h_mCWR_SE[KEY_SE]->GetNbinsX() << ", bin width = " << h_mCWR_SE[KEY_SE]->GetBinWidth(1) << ", low edge = " << h_mCWR_SE[KEY_SE]->GetBinLowEdge(1) << ", high edge = " << h_mCWR_SE[KEY_SE]->GetBinLowEdge(h_mCWR_SE[KEY_SE]->GetNbinsX()) + h_mCWR_SE[KEY_SE]->GetBinWidth(1) << endl;        
      //cout << "nbins GAVIN = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetNbinsX() << ", bin width = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinWidth(1) << ", low edge = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinLowEdge(1) << ", high edge = " << h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinLowEdge(h_mGAVIN_SE[KEY_GAVIN_SE]->GetNbinsX()) + h_mGAVIN_SE[KEY_GAVIN_SE]->GetBinWidth(1) << endl;        
    }
    string KEY_ratio_SE_Int = Form("Centrality_9_CosThetaStar_%d_SE_Ratio",i_theta);
    string KEY_ratio_ME_Int = Form("Centrality_9_CosThetaStar_%d_ME_Ratio",i_theta);
    h_mRatio_SE[KEY_ratio_SE_Int] = (TH1F*) h_mCWR_SE[KEY_SE_Int]->Clone(KEY_ratio_SE_Int.c_str());
    h_mRatio_SE[KEY_ratio_SE_Int]->Divide(h_mGAVIN_SE[KEY_GAVIN_SE_Int]);
  }

  TCanvas *c_SE = new TCanvas("c_SE","c_SE",10,10,1200,800);
  c_SE->Divide(3,2);
  for(int i = 0; i < pt_rebin_local; i++)
  {
    c_SE->cd(i+1)->SetLeftMargin(0.15);
    c_SE->cd(i+1)->SetBottomMargin(0.15);
    c_SE->cd(i+1)->SetTicks(1,1);
    c_SE->cd(i+1)->SetGrid(0,0);
  }

  string outputname = Form("figures/SE_GAVIN_CWR_COMPARISON_27GEV_RUN18.pdf");
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop = Form("%s]",outputname.c_str());

  c_SE->Print(outputstart.c_str());

  for(int i_theta = 0; i_theta < 7; i_theta++)
  {
    for(int pt_bin = 1; pt_bin < pt_rebin_local; pt_bin++)
    {
      c_SE->cd(pt_bin);
      string KEY_CWR_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_CWR",pt_bin,i_theta);
      //string KEY_CWR_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_CWR",pt_bin,i_theta);
      string KEY_GAVIN_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",pt_bin,i_theta);
      //string KEY_GAVIN_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",pt_bin,i_theta);
      h_mGAVIN_SE[KEY_GAVIN_SE]->SetTitle(Form("27GeV R18, 20-60, %1.1f<p_{T}<%1.1f, %d/7<|cos(#theta*)|<%d/7",pt_low_local[pt_bin],pt_up_local[pt_bin],i_theta,i_theta+1));
      h_mGAVIN_SE[KEY_GAVIN_SE]->GetYaxis()->SetTitle("Counts");
      h_mGAVIN_SE[KEY_GAVIN_SE]->GetXaxis()->SetTitle("M(K^{+},K^{-}) GeV");
      h_mCWR_SE[KEY_CWR_SE]->SetMarkerStyle(20);
      h_mCWR_SE[KEY_CWR_SE]->SetMarkerColor(kOrange+7);
      h_mCWR_SE[KEY_CWR_SE]->SetLineColor(kOrange+7);

      h_mGAVIN_SE[KEY_GAVIN_SE]->SetMarkerStyle(20);
      h_mGAVIN_SE[KEY_GAVIN_SE]->SetMarkerColor(kBlue);
      h_mGAVIN_SE[KEY_GAVIN_SE]->SetLineColor(kBlue);

      h_mGAVIN_SE[KEY_GAVIN_SE]->Draw("pE");
      h_mCWR_SE[KEY_CWR_SE]->Draw("pE same");

      TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
      leg1->AddEntry(h_mCWR_SE[KEY_CWR_SE],"SE CWR","p");
      leg1->AddEntry(h_mGAVIN_SE[KEY_GAVIN_SE],"SE GAVIN Scaled","p");
      leg1->Draw("same");
    }
    c_SE->cd(6);
    string KEY_CWR_SE = Form("Centrality_9_CosThetaStar_%d_SE_Rebinned_CWR",i_theta);
    string KEY_GAVIN_SE = Form("Centrality_9_CosThetaStar_%d_SE_Rebinned_GAVIN",i_theta);
    h_mGAVIN_SE[KEY_GAVIN_SE]->SetTitle(Form("27GeV R18, 20-60, 1.2<p_{T}<5.4, %d/7<|cos(#theta*)|<%d/7",i_theta,i_theta+1));
    h_mGAVIN_SE[KEY_GAVIN_SE]->GetYaxis()->SetTitle("Counts");
    h_mGAVIN_SE[KEY_GAVIN_SE]->GetXaxis()->SetTitle("M(K^{+},K^{-}) GeV");
    h_mCWR_SE[KEY_CWR_SE]->SetMarkerStyle(20);
    h_mCWR_SE[KEY_CWR_SE]->SetMarkerColor(kOrange+7);
    h_mCWR_SE[KEY_CWR_SE]->SetLineColor(kOrange+7);

    h_mGAVIN_SE[KEY_GAVIN_SE]->SetMarkerStyle(20);
    h_mGAVIN_SE[KEY_GAVIN_SE]->SetMarkerColor(kBlue);
    h_mGAVIN_SE[KEY_GAVIN_SE]->SetLineColor(kBlue);

    h_mGAVIN_SE[KEY_GAVIN_SE]->Draw("pE");
    h_mCWR_SE[KEY_CWR_SE]->Draw("pE same");

    TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
    leg1->AddEntry(h_mCWR_SE[KEY_CWR_SE],"SE CWR","p");
    leg1->AddEntry(h_mGAVIN_SE[KEY_GAVIN_SE],"SE GAVIN Scaled","p");
    leg1->Draw("same");
    
    c_SE->Update();
    c_SE->Print(outputname.c_str());
  }
  c_SE->Print(outputstop.c_str());


  outputname = Form("figures/SERATIO_GAVIN_CWR_COMPARISON_27GEV_RUN18.pdf");
  outputstart = Form("%s[",outputname.c_str());
  outputstop = Form("%s]",outputname.c_str());

  c_SE->Print(outputstart.c_str());

  for(int i_theta = 0; i_theta < 7; i_theta++)
  {
    for(int pt_bin = 1; pt_bin < pt_rebin_local; pt_bin++)
    {
      c_SE->cd(pt_bin);
      string KEY_ratio_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_SE_Ratio",pt_bin,i_theta);
      h_mRatio_SE[KEY_ratio_SE]->SetTitle(Form("27GeV R18, 20-60, %1.1f<p_{T}<%1.1f, %d/7<|cos(#theta*)|<%d/7",pt_low_local[pt_bin],pt_up_local[pt_bin],i_theta,i_theta+1));
      double min = h_mRatio_SE[KEY_ratio_SE]->GetMinimum(0.01);
      double max = h_mRatio_SE[KEY_ratio_SE]->GetMaximum();
      h_mRatio_SE[KEY_ratio_SE]->GetYaxis()->SetRangeUser(min*0.975,max*1.025);
      h_mRatio_SE[KEY_ratio_SE]->GetYaxis()->SetTitle("CWR/(GAVIN scaled)");
      h_mRatio_SE[KEY_ratio_SE]->GetXaxis()->SetTitle("M(K^{+},K^{-}) GeV");
      h_mRatio_SE[KEY_ratio_SE]->SetMarkerStyle(20);
      h_mRatio_SE[KEY_ratio_SE]->SetMarkerColor(kOrange+7);
      h_mRatio_SE[KEY_ratio_SE]->SetLineColor(kOrange+7);

      h_mRatio_SE[KEY_ratio_SE]->Draw("pE");
    }
    c_SE->cd(6);
    string KEY_ratio_SE = Form("Centrality_9_CosThetaStar_%d_SE_Ratio",i_theta);
    h_mRatio_SE[KEY_ratio_SE]->SetTitle(Form("27GeV R18, 20-60, 1.2<p_{T}<5.4, %d/7<|cos(#theta*)|<%d/7",i_theta,i_theta+1));
    double min = h_mRatio_SE[KEY_ratio_SE]->GetMinimum(0.01);
    double max = h_mRatio_SE[KEY_ratio_SE]->GetMaximum();
    h_mRatio_SE[KEY_ratio_SE]->GetYaxis()->SetRangeUser(min*0.975,max*1.025);
    h_mRatio_SE[KEY_ratio_SE]->GetYaxis()->SetTitle("CWR/(GAVIN scaled)");
    h_mRatio_SE[KEY_ratio_SE]->GetXaxis()->SetTitle("M(K^{+},K^{-}) GeV");
    h_mRatio_SE[KEY_ratio_SE]->SetMarkerStyle(20);
    h_mRatio_SE[KEY_ratio_SE]->SetMarkerColor(kOrange+7);
    h_mRatio_SE[KEY_ratio_SE]->SetLineColor(kOrange+7);
    h_mRatio_SE[KEY_ratio_SE]->Draw("pE");

    c_SE->Update();
    c_SE->Print(outputname.c_str());
  }
  c_SE->Print(outputstop.c_str());
 
  string outputfile = Form("../output/AuAu%s/%s/CWRPrepped.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
  {
    string KEY_SE_Int = Form("Centrality_9_CosThetaStar_%d_SE_Rebinned_CWR",i_theta);
    string KEY_ME_Int = Form("Centrality_9_CosThetaStar_%d_ME_Rebinned_CWR",i_theta);
    h_mCWR_SE[KEY_SE_Int]->Write();
    h_mCWR_ME[KEY_ME_Int]->Write();
    
  }
  File_OutPut->Close();
}
