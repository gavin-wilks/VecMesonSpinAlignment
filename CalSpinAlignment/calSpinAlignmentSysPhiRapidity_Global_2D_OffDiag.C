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

#ifndef _PlotQA_
#define _PlotQA_  1
#endif


void calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag(int energy = 4, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1", int frameopt = 0)
{
  std::string frame = "Global";
  if(frameopt == 1) frame = "Helicity";

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000); 
  
  std::string EP[2] = {"1st","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/%s2DRapidity_InvMassSubBg_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/%s2DRapidity_InvMassSubBg_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  //File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
 
  //TGraMap g_mChiNDF;
  //TGraMap g_mPValue;
  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  for(int  i_eta = 0; i_eta < vmsa::eta_total; i_eta++) // pt loop
  {
    for(int i_pt = 0; i_pt < vmsa::pt_rebin_y_2D; i_pt++)
    {
      for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++) // Centrality loop
      {
        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
          {
            if( i_dca != 0 && i_sig != 0 ) continue;
            for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
            {
              string KEY_InteTheta = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              for(int i_phi = vmsa::CTS_start; i_phi < 12/*2*vmsa::CTS_stop*/; i_phi++) // cos(theta*) loop
              {
                for(int i_thetah = vmsa::CTS_start; i_thetah < 9/*2*vmsa::CTS_stop*/; i_thetah++) // cos(theta*) loop
                {
                  string KEY = Form("pt_%d_eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str());

                  if(i_thetah == vmsa::CTS_start && i_phi == vmsa::CTS_start) h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)h_mMass[KEY]->Clone(KEY_InteTheta.c_str());
                  else h_mMass_InteTheta[KEY_InteTheta]->Add(h_mMass[KEY],1.0);
                }
              }
              for(int i_poly = 0; i_poly < 3; i_poly++)
              {
                string KEY_InteTheta_Poly = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                //string KEY_Poly = Form("Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                //if(i_eta == vmsa::pt_rebin_first[energy]) g_mChiNDF[KEY_Poly] = new TGraphAsymmErrors();             
                //if(i_eta == vmsa::pt_rebin_first[energy]) g_mPValue[KEY_Poly] = new TGraphAsymmErrors();             

                TF1 *f_bw; 
                if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
                if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
                if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
                f_bw->SetParameter(0,vmsa::InvMass[pid]);
                f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.005,vmsa::InvMass[pid]+0.005);
                f_bw->SetParameter(1,vmsa::Width[pid]);
                f_bw->SetParameter(2,1.0);
                float norm = h_mMass_InteTheta[KEY_InteTheta]->GetMaximum()/f_bw->GetMaximum();
                f_bw->SetParameter(2,norm);
                f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
                h_mMass_InteTheta[KEY_InteTheta]->Fit(f_bw,"MQNR");
                Par_InteTheta[KEY_InteTheta_Poly].clear();
                Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(0)));
                Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(1)));
                Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(2)));
                Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(3)));
                Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(4)));
                if(i_poly >= 1) Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(5)));
                if(i_poly >= 2) Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(6)));
              }
            }
          }
        }
      }
    }
  }

  // extract counts vs. pT with diffenretial integration ranges and methods
  vecFMap Par;
  TH2FMap h_mCounts;
  TH2FMap h_mRatio;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;
  TGraMap g_mReal;
  TGraMap g_mImag;
  TGraMap g_mReRho1n1;
  TGraMap g_mImRho1n1;

  bool skip[10][10][vmsa::eta_total][9][12] = {false};
  for(int i_pt = 0; i_pt < vmsa::pt_rebin_y_2D; i_pt++)
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;
          for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
          {
            for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
            {
              for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
              {
                for(int i_poly = 0; i_poly < 3; i_poly++)
                {
                  string KEY_rho = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  string KEY_real = Form("realRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  string KEY_imag = Form("imagRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  string KEY_rerho1n1 = Form("rerho1n1Raw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  string KEY_imrho1n1 = Form("imrho1n1Raw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  //cout << KEY_rho << endl;
                  g_mRho[KEY_rho] = new TGraphAsymmErrors();
                  g_mReal[KEY_real] = new TGraphAsymmErrors();
                  g_mImag[KEY_imag] = new TGraphAsymmErrors();
                  g_mReRho1n1[KEY_rerho1n1] = new TGraphAsymmErrors();
                  g_mImRho1n1[KEY_imrho1n1] = new TGraphAsymmErrors();
                  for(int  i_eta = 0; i_eta < vmsa::eta_total; ++i_eta) // pt loop
                  {
                    string KEY_counts = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                    h_mCounts[KEY_counts] = new TH2F(KEY_counts.c_str(),KEY_counts.c_str(),9,-1.0,1.0,12,0.0,2.0*TMath::Pi());

                    string KEY_InteTheta = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                    string KEY_InteTheta_Poly = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                    //g_mChiNDF[KEY_InteTheta_Poly] = new TGraphAsymmErrors();             
                    //g_mPValue[KEY_InteTheta_Poly] = new TGraphAsymmErrors();             
                    for(int i_phi = vmsa::CTS_start; i_phi < 12/*2*vmsa::CTS_stop*/; ++i_phi) // cos(theta*) loop
                    {
                      cout << "i_phi = " << i_phi << endl;
                      for(int i_thetah = vmsa::CTS_start; i_thetah < 9/*2*vmsa::CTS_stop*/; ++i_thetah) // cos(theta*) loop
                      {
                        string KEY = Form("pt_%d_eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                        string KEY_Poly = Form("pt_%d_eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_eta,i_cent,i_thetah,i_phi,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                        TF1 *f_bw;
                        if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
                        if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
                        if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
                        f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_Poly][0]);
                        f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.005,vmsa::InvMass[pid]+0.005);
                        f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_Poly][1]);
                        f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_Poly][2]/108.0);
                        f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_Poly][3]/108.0);
                        f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_Poly][4]/108.0);
                        if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_Poly][5]);
                        if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_Poly][6]);
                        //f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta][5]);
                        f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
                        //if(h_mMass[KEY]->GetEntries() <= 1) continue;
                    
                        cout << KEY << endl;
                        cout << "Number of entries = " << h_mMass[KEY]->GetEntries() << endl;
                        TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"MQNRS");
                        if(!result.Get()) 
                        {
                          skip[i_pt][i_cent][i_eta][i_thetah][i_phi] = true;
                          cout << "WE SKIPPED THIS" << endl;
                          continue;
                        } 
                        TF1 *f_bg;
                        if(i_poly == 0) f_bg = new TF1("f_bg",Poly, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
                        if(i_poly == 1) f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
                        if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
                        f_bg->SetParameter(0,f_bw->GetParameter(3));
                        f_bg->SetParameter(1,f_bw->GetParameter(4));
                        if(i_poly >= 1) f_bg->SetParameter(2,f_bw->GetParameter(5));
                        if(i_poly >= 2) f_bg->SetParameter(3,f_bw->GetParameter(6));
                        f_bg->SetParError(0,f_bw->GetParError(3));
                        f_bg->SetParError(1,f_bw->GetParError(4));
                        if(i_poly >= 1) f_bg->SetParError(2,f_bw->GetParError(5));
                        if(i_poly >= 2) f_bg->SetParError(3,f_bw->GetParError(6));

                        double params1[2] = {0.0};//,result->GetParams()[5];
                        if(i_poly == 0) 
                        {
                          params1[0] = result->GetParams()[3];
                          params1[1] = result->GetParams()[4];
                        }
                        double params2[3] = {0.0};
                        if(i_poly == 1) 
                        {
                          params2[0] = result->GetParams()[3];
                          params2[1] = result->GetParams()[4];
                          params2[2] = result->GetParams()[5];
                        }
                        double params3[4] = {0.0};
                        if(i_poly == 2) 
                        {
                          params3[0] = result->GetParams()[3];
                          params3[1] = result->GetParams()[4];
                          params3[2] = result->GetParams()[5];
                          params3[2] = result->GetParams()[6];
                        }

                        TMatrixDSym covArr1(2);
                        TMatrixDSym covArr2(3);
                        TMatrixDSym covArr3(4);
                        if(i_poly == 0)
                        {
                          covArr1(0,0) = result->GetCovarianceMatrix()(3,3);
                          covArr1(0,1) = result->GetCovarianceMatrix()(3,4);
                          covArr1(1,0) = result->GetCovarianceMatrix()(4,3);
                          covArr1(1,1) = result->GetCovarianceMatrix()(4,4);
                        }
                        if(i_poly == 1)
                        {
                          covArr2(0,0) = result->GetCovarianceMatrix()(3,3);
                          covArr2(0,1) = result->GetCovarianceMatrix()(3,4);
                          covArr2(0,2) = result->GetCovarianceMatrix()(3,5);
                          covArr2(1,0) = result->GetCovarianceMatrix()(4,3);
                          covArr2(1,1) = result->GetCovarianceMatrix()(4,4);
                          covArr2(1,2) = result->GetCovarianceMatrix()(4,5);
                          covArr2(2,0) = result->GetCovarianceMatrix()(5,3);
                          covArr2(2,1) = result->GetCovarianceMatrix()(5,4);
                          covArr2(2,2) = result->GetCovarianceMatrix()(5,5);
                        }
                        if(i_poly == 2)
                        {
                          covArr3(0,0) = result->GetCovarianceMatrix()(3,3);
                          covArr3(0,1) = result->GetCovarianceMatrix()(3,4);
                          covArr3(0,2) = result->GetCovarianceMatrix()(3,5);
                          covArr3(0,3) = result->GetCovarianceMatrix()(3,6);
                          covArr3(1,0) = result->GetCovarianceMatrix()(4,3);
                          covArr3(1,1) = result->GetCovarianceMatrix()(4,4);
                          covArr3(1,2) = result->GetCovarianceMatrix()(4,5);
                          covArr3(1,3) = result->GetCovarianceMatrix()(4,6);
                          covArr3(2,0) = result->GetCovarianceMatrix()(5,3);
                          covArr3(2,1) = result->GetCovarianceMatrix()(5,4);
                          covArr3(2,2) = result->GetCovarianceMatrix()(5,5);
                          covArr3(2,3) = result->GetCovarianceMatrix()(5,6);
                          covArr3(3,0) = result->GetCovarianceMatrix()(6,3);
                          covArr3(3,1) = result->GetCovarianceMatrix()(6,4);
                          covArr3(3,2) = result->GetCovarianceMatrix()(6,5);
                          covArr3(3,3) = result->GetCovarianceMatrix()(6,6);
                        }

                        float bin_width = h_mMass[KEY]->GetBinWidth(1);
                        float Inte_start = Par_InteTheta[KEY_InteTheta_Poly][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]-0.5*bin_width;
                        float Inte_stop  = Par_InteTheta[KEY_InteTheta_Poly][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]+0.5*bin_width;
                        float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
                        float errors_bg; 
                        if(i_poly == 0) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params1,covArr1.GetMatrixArray())/bin_width;
                        if(i_poly == 1) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params2,covArr2.GetMatrixArray())/bin_width;
                        if(i_poly == 2) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params3,covArr3.GetMatrixArray())/bin_width;

                        //float cos_mean = 1/14.0+float(i_theta-7)/7.0;
                        //float chi2NDF = float(f_bw->GetChisquare())/float(f_bw->GetNDF()); 
                        //float pvalue = TMath::Prob(f_bw->GetChisquare(),f_bw->GetNDF());
                        ////cout << "Poly " << i_poly + 1 << ", cos = " << cos_mean << ", chi2 = " << f_bw->GetChisquare() << ", NDF = " << f_bw->GetNDF() << ", pvalue = " << pvalue << endl;
                        //g_mChiNDF[KEY_InteTheta_Poly]->SetPoint(i_theta,cos_mean,chi2NDF);
                        //g_mPValue[KEY_InteTheta_Poly]->SetPoint(i_theta,cos_mean,pvalue);

                        float phi_bin_center = TMath::Pi()*float(i_phi+0.5)/6.0;
                        float bin_centerh = float(i_thetah-4)/4.5;
                        if(i_method == 0)
                        {
                          int bin_start = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta_Poly][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]);
                          int bin_stop  = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta_Poly][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]);
                          float counts = 0.0;
                          float errors = 0.0;
                          for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
                          {
                            counts += h_mMass[KEY]->GetBinContent(i_bin);
                            errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
                          }
                          h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_centerh,phi_bin_center),counts-counts_bg);
                          h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_centerh,phi_bin_center),TMath::Sqrt(errors+errors_bg*errors_bg));
                        }
                        if(i_method == 1)
                        {
                          //float bin_width = h_mMass[KEY]->GetBinWidth(1);
                          //float Inte_start = Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
                          //float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
                          float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
                          float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
                          h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_centerh,phi_bin_center),counts_bw-counts_bg);
                          h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_centerh,phi_bin_center),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));
                        }
                        Par[KEY_Poly].clear();
                        Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(0)));
                        Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(1)));
                        Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(2)));
                        Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(3)));
                        Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(4)));
                        if(i_poly >= 1) Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(5)));
                        if(i_poly >= 2) Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(6)));
                      }
                    }
                    //float pt_mean = float(float(i_eta)-4.5)/5.0;
                    float pt_mean = (vmsa::eta_bin[i_eta]+vmsa::eta_bin[i_eta+1])/2.0;

                    TF2 *f_rho = new TF2("f_rho",SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
                    f_rho->SetParameter(0,0.33);
                    f_rho->SetParameter(1,0.0);
                    f_rho->SetParameter(2,0.0);
                    f_rho->SetParameter(3,0.0);
                    f_rho->SetParameter(4,0.0);
                    f_rho->SetParameter(5,h_mCounts[KEY_counts]->GetMaximum());
                    h_mCounts[KEY_counts]->Fit(f_rho,"NMRI");
                    Par_rhoFit[KEY_counts].clear();
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(2)));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(3)));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(4)));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(5)));
                    g_mRho[KEY_rho]->SetPoint(i_eta,pt_mean,f_rho->GetParameter(0)-1./3.);
                    g_mRho[KEY_rho]->SetPointError(i_eta,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
                    g_mReal[KEY_real]->SetPoint(i_eta,pt_mean,f_rho->GetParameter(1));
                    g_mReal[KEY_real]->SetPointError(i_eta,0.0,0.0,f_rho->GetParError(1),f_rho->GetParError(1));
                    g_mImag[KEY_imag]->SetPoint(i_eta,pt_mean,f_rho->GetParameter(2));
                    g_mImag[KEY_imag]->SetPointError(i_eta,0.0,0.0,f_rho->GetParError(2),f_rho->GetParError(2));
                    g_mReRho1n1[KEY_rerho1n1]->SetPoint(i_eta,pt_mean,f_rho->GetParameter(3));
                    g_mReRho1n1[KEY_rerho1n1]->SetPointError(i_eta,0.0,0.0,f_rho->GetParError(3),f_rho->GetParError(3));
                    g_mImRho1n1[KEY_imrho1n1]->SetPoint(i_eta,pt_mean,f_rho->GetParameter(4));
                    g_mImRho1n1[KEY_imrho1n1]->SetPointError(i_eta,0.0,0.0,f_rho->GetParError(4),f_rho->GetParError(4));
                    //cout << "for rho00 calculation, pt_mean = " << pt_mean << ", rho = " << f_rho->GetParameter(0) << " +/- " << f_rho->GetParError(0) << endl;
                    
                  }
                }
              }
            }
          }
        }
      }
  }
  
#if _PlotQA_
  for(int i_poly = 0; i_poly < 1/*3*/; i_poly++)
  {
    for(int i_norm = 0; i_norm < 1/*3*/; i_norm++)
    {
      for(int i_pt = 0; i_pt < vmsa::pt_rebin_y_2D; i_pt++)
      {
        for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; ++i_cent) // Centrality loop
        {
          string outputname = Form("./figures/%s/%s/rapiditystudy/Global2D_allThetaYields%s_%s_Order%d_cent%d_%s_pTdependence_Norm%d_Poly%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_cent,etamode.c_str(),i_norm,i_poly+1);
          string output_start = Form("%s[",outputname.c_str());
          
          TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
          TCanvas *c_diff2d = new TCanvas("c_diff2d","c_diff2d",10,10,500,500);
          c_diff->Print(output_start.c_str());
          for(int  i_eta = 0; i_eta < vmsa::eta_total; ++i_eta) // pt loop
          {
            for(int i_phi = 0; i_phi < 12; i_phi++)
            { 
              c_diff->Clear();
              c_diff->Divide(4,3);
              for(int i_thetah = 0; i_thetah < 9; ++i_thetah)
              {
                cout << "i_cent = " << i_cent << ", i_eta = " << i_eta << ", i_thetah = " << i_thetah << ", i_phi = " << i_phi << ", skip? " << skip[i_cent][i_eta][i_thetah][i_phi] << endl;
                if(skip[i_cent][i_eta][i_thetah][i_phi]) continue;
                c_diff->cd(i_thetah+1);
                c_diff->cd(i_thetah+1)->SetLeftMargin(0.15);
                c_diff->cd(i_thetah+1)->SetBottomMargin(0.15);
                c_diff->cd(i_thetah+1)->SetTicks(1,1);
                c_diff->cd(i_thetah+1)->SetGrid(0,0);
                if(i_thetah < 10)
                {
                  string KEY_QA = Form("eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_cent,i_thetah,i_phi,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm);
                  //cout << KEY_QA << endl;
                  h_mMass[KEY_QA]->Print();
                  h_mMass[KEY_QA]->SetTitle("");
                  h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
                  h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
                  h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
                  h_mMass[KEY_QA]->SetTitle(Form("%.2f<y<%.2f phiH=%d, CTH=%d",vmsa::eta_rebinbin[i_eta],vmsa::eta_rebinbin[i_eta+1],i_phi,i_thetah));
                  h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
                  h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
                  h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
  
                  h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
                  h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
                  h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
                  h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
                  h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
                  h_mMass[KEY_QA]->GetYaxis()->CenterTitle();
  
                  h_mMass[KEY_QA]->SetMarkerStyle(24);
                  h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
                  h_mMass[KEY_QA]->SetMarkerSize(1.2);
                  h_mMass[KEY_QA]->Draw("pE");
                  PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  
                  //string KEY_InteTheta = Form("eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,1/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
                  //float x1_0 = Par_InteTheta[KEY_InteTheta][0] - 2.0*Par_InteTheta[KEY_InteTheta][1];
                  //float x2_0 = Par_InteTheta[KEY_InteTheta][0] + 2.0*Par_InteTheta[KEY_InteTheta][1];
                  //float x1_1 = Par_InteTheta[KEY_InteTheta][0] - 2.5*Par_InteTheta[KEY_InteTheta][1];
                  //float x2_1 = Par_InteTheta[KEY_InteTheta][0] + 2.5*Par_InteTheta[KEY_InteTheta][1];
                  //float x1_2 = Par_InteTheta[KEY_InteTheta][0] - 3.0*Par_InteTheta[KEY_InteTheta][1];
                  //float x2_2 = Par_InteTheta[KEY_InteTheta][0] + 3.0*Par_InteTheta[KEY_InteTheta][1];
                  //float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(Par_InteTheta[KEY_InteTheta][0]));
                  //h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(-0.3*y,1.2*y);
                  //float ymin = h_mMass[KEY_QA]->GetMinimum();
                  //PlotLine(x1_0,x1_0,-0.3*y,y,4,2,2);
                  //PlotLine(x2_0,x2_0,-0.3*y,y,4,2,2);
                  //PlotLine(x1_1,x1_1,-0.3*y,y,1,2,2);
                  //PlotLine(x2_1,x2_1,-0.3*y,y,1,2,2);
                  //PlotLine(x1_2,x1_2,-0.3*y,y,4,2,2);
                  //PlotLine(x2_2,x2_2,-0.3*y,y,4,2,2);
                }
                if(i_thetah == 10)
                {
                  string KEY_InteTheta_QA = Form("eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm);
                  string KEY_InteTheta_QA_Poly = Form("eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_eta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->Print();
                  h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
                  h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle(Form("%.2f<y<%.2f, Full phiH, Full CTG",vmsa::eta_rebinbin[i_eta],vmsa::eta_rebinbin[i_eta+1]));
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
  
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Yields");
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();
  
                  h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
                  h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
                  PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
                  TF1 *f_bw;
                  cout << KEY_InteTheta_QA_Poly << endl;
                  if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
                  if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
                  if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
                  f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA_Poly][0]);
                  f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA_Poly][1]);
                  f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA_Poly][2]);
                  f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA_Poly][3]);
                  f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA_Poly][4]);
                  if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA_Poly][5]);
                  if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_QA_Poly][6]);
                  f_bw->SetLineColor(2);
                  f_bw->SetLineStyle(2);
                  f_bw->SetLineWidth(2);
                  f_bw->Draw("l same");
                }
              }
  
              for(int i_thetah = vmsa::CTS_start; i_thetah < 10/*&2*vmsa::CTS_stop*/; ++i_thetah)
              {
                if(skip[i_cent][i_eta][i_thetah][i_phi]) continue;
                c_diff->cd(i_thetah+1);
                string KEY_QA = Form("eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_eta,i_cent,i_thetah,i_phi,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm);
                string KEY_QA_Poly = Form("eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_eta,i_cent,i_thetah,i_phi,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                cout << KEY_QA << endl;
                cout << KEY_QA_Poly << endl;
                TF1 *f_bw;
                if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
                if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
                if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
                f_bw->SetParameter(0,Par[KEY_QA_Poly][0]);
                f_bw->SetParameter(1,Par[KEY_QA_Poly][1]);
                f_bw->SetParameter(2,Par[KEY_QA_Poly][2]);
                f_bw->SetParameter(3,Par[KEY_QA_Poly][3]);
                f_bw->SetParameter(4,Par[KEY_QA_Poly][4]);
                if(i_poly >= 1) f_bw->SetParameter(5,Par[KEY_QA_Poly][5]);
                if(i_poly >= 2) f_bw->SetParameter(6,Par[KEY_QA_Poly][6]);
                f_bw->SetLineColor(kOrange+7);
                f_bw->SetLineStyle(1);
                f_bw->SetLineWidth(2);
                f_bw->Draw("l same");
  
                TF1 *f_bg;
                if(i_poly == 0) f_bg = new TF1("f_bg",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
                if(i_poly == 1) f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
                if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
                f_bg->SetParameter(0,Par[KEY_QA_Poly][3]);
                f_bg->SetParameter(1,Par[KEY_QA_Poly][4]);
                if(i_poly >= 1) f_bg->SetParameter(2,Par[KEY_QA_Poly][5]);
                if(i_poly >= 2) f_bg->SetParameter(3,Par[KEY_QA_Poly][6]);
                f_bg->SetLineColor(kBlue);
                f_bg->SetLineStyle(2);
                f_bg->SetLineWidth(2);
                f_bg->Draw("l same");
  
                TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
                leg1->AddEntry(h_mMass[KEY_QA],"data","p");
                leg1->AddEntry(f_bw,"sig+res","l");
                leg1->AddEntry(f_bg,"res","l");
                leg1->Draw("same");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
            }
            c_diff2d->cd();
            c_diff2d->cd()->SetLeftMargin(0.15);
            c_diff2d->cd()->SetRightMargin(0.15);
            c_diff2d->cd()->SetBottomMargin(0.15);
            c_diff2d->cd()->SetTicks(1,1);
            c_diff2d->cd()->SetGrid(0,0);
            for(int i_sigma = 0/*vmsa::Sig_start*/; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
            {
              for(int i_method = 1/*vmsa::Method_start*/; i_method < 2/*vmsa::Method_stop*/; ++i_method)
              {
                string KEY_counts_QA = Form("eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_eta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                cout << KEY_counts_QA << endl;
                h_mCounts[KEY_counts_QA]->SetStats(0);
                h_mCounts[KEY_counts_QA]->Print();
                h_mCounts[KEY_counts_QA]->SetTitle(Form("Raw Yields, %.1f<y<%.1f",vmsa::eta_rebinbin[i_eta],vmsa::eta_rebinbin[i_eta+1]));
                h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
                h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
                h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta)");
                h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
                h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
                h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();
  
                h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
                h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
                h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("#beta");
                h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
                h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
                h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();
  
                h_mCounts[KEY_counts_QA]->Draw("colz");
              }
            }
            c_diff2d->Update();
            c_diff2d->Print(outputname.c_str());
          }
          string output_stop = Form("%s]",outputname.c_str());
          c_diff->Print(output_stop.c_str()); // close pdf file
        }
      }
    }
  }
#endif

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.2,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);


  TCanvas *c_spin = new TCanvas("c_spin","c_spin",10,10,900,600);
  c_spin->Divide(3,2);
  for(int i = 0; i < 6; i++)
  {
    c_spin->cd(i+1);
    c_spin->cd(i+1)->SetLeftMargin(0.15);
    c_spin->cd(i+1)->SetBottomMargin(0.15);
    c_spin->cd(i+1)->SetTicks(1,1);
    c_spin->cd(i+1)->SetGrid(0,0);
  }
   
  string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),          0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_real = Form("realRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),        0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_imag = Form("imagRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),        0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_rerho1n1 = Form("rerho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_imrho1n1 = Form("imrho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

  c_spin->cd(1);
  g_mRho[KEY_rho]->GetYaxis()->SetTitle("#rho_{00}");
  g_mRho[KEY_rho]->GetXaxis()->SetTitle("y");
  g_mRho[KEY_rho]->SetMarkerStyle(20); 
  g_mRho[KEY_rho]->Draw("APE");

  c_spin->cd(2);
  g_mReal[KEY_real]->GetYaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
  g_mReal[KEY_real]->GetXaxis()->SetTitle("y");
  g_mReal[KEY_real]->SetMarkerStyle(20); 
  g_mReal[KEY_real]->Draw("APE");

  c_spin->cd(3);
  g_mImag[KEY_imag]->GetYaxis()->SetTitle("Im(#rho_{10})-Im(#rho_{0-1})");
  g_mImag[KEY_imag]->GetXaxis()->SetTitle("y");
  g_mImag[KEY_imag]->SetMarkerStyle(20); 
  g_mImag[KEY_imag]->Draw("APE");

  c_spin->cd(4);
  g_mReRho1n1[KEY_rerho1n1]->GetYaxis()->SetTitle("Re(#rho_{1-1})");
  g_mReRho1n1[KEY_rerho1n1]->GetXaxis()->SetTitle("y");
  g_mReRho1n1[KEY_rerho1n1]->SetMarkerStyle(20); 
  g_mReRho1n1[KEY_rerho1n1]->Draw("APE");

  c_spin->cd(5);
  g_mImRho1n1[KEY_imrho1n1]->GetYaxis()->SetTitle("Im(#rho_{1-1})");
  g_mImRho1n1[KEY_imrho1n1]->GetXaxis()->SetTitle("y");
  g_mImRho1n1[KEY_imrho1n1]->SetMarkerStyle(20); 
  g_mImRho1n1[KEY_imrho1n1]->Draw("APE");

  string outputname = Form("./figures/%s/%s/rapiditystudy/Global2D_ExtractedSpinDensityElements_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  c_spin->Print(outputname.c_str());
  

  //string outputfile = Form("../output/AuAu%s/%s/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string outputfile = Form("../output/AuAu%s/%s/Global2DRapidity_RawPhiYSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/Global2DRapidity_RawPhiYSys_%s_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
	  {
	    for(int i_method = 1/*vmsa::Method_start*/; i_method < vmsa::Method_stop; ++i_method)
	    {
              for(int i_poly = 0; i_poly < 1/*3*/; i_poly++)
              {
	        string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	        string KEY_real = Form("realRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	        string KEY_imag = Form("imagRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	        string KEY_rerho1n1 = Form("rerho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	        string KEY_imrho1n1 = Form("imrho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	        g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
	        g_mRho[KEY_rho]->Write();
	        g_mReal[KEY_real]->SetName(KEY_real.c_str());
	        g_mReal[KEY_real]->Write();
	        g_mImag[KEY_imag]->SetName(KEY_imag.c_str());
	        g_mImag[KEY_imag]->Write();
	        g_mReRho1n1[KEY_rerho1n1]->SetName(KEY_rerho1n1.c_str());
	        g_mReRho1n1[KEY_rerho1n1]->Write();
	        g_mImRho1n1[KEY_imrho1n1]->SetName(KEY_imrho1n1.c_str());
	        g_mImRho1n1[KEY_imrho1n1]->Write();
	        for(int i_eta = 0; i_eta < vmsa::eta_rebintotal; ++i_eta) // pt loop
	        {
	          string KEY_counts = Form("eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	          h_mCounts[KEY_counts]->Print();
	          h_mCounts[KEY_counts]->Write();
	        }
              }
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}
