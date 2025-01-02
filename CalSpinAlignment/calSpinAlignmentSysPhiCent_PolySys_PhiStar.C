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


void calSpinAlignmentSysPhiCent_PolySys_PhiStar(int energy = 4, int pid = 0, int year = 0, bool random3D = false, int centQA = 4, int ptQA = 0, int dcaQA = 0, int nsigQA = 0, int order = 2, std::string etamode = "eta1_eta1")
{
  std::string EP[2] = {"","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_InvMassSubBgCent_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_InvMassSubBgCent_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBgCent.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;

  TGraMap g_mChiNDF;
  TGraMap g_mPValue;

  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter

  for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;
          for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
          {
            string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
            for(int i_theta = 0; i_theta < 10; i_theta++) // cos(theta*) loop
            {
              string KEY = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str())->Clone();
              h_mMass[KEY]->Rebin(2);

              if(i_theta == 0) h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)h_mMass[KEY]->Clone(KEY_InteTheta.c_str());
              else h_mMass_InteTheta[KEY_InteTheta]->Add(h_mMass[KEY],1.0);
              //cout << KEY << endl;
            }
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_InteTheta_Poly = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
              TF1 *f_bw;
              if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
              if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
              if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
              f_bw->SetParameter(0,vmsa::InvMass[pid]);
              f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.003,vmsa::InvMass[pid]+0.003);
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
  

//#if _PlotQA_
//  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,900);
//  c_diff->Divide(3,3);
//  for(int i_theta = 0; i_theta < 9; ++i_theta)
//  {
//    c_diff->cd(i_theta+1);
//    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
//    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
//    c_diff->cd(i_theta+1)->SetTicks(1,1);
//    c_diff->cd(i_theta+1)->SetGrid(0,0);
//    if(i_theta < 10)
//    {
//      string KEY_QA = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",ptQA,centQA,i_theta,dcaQA,nsigQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//      h_mMass[KEY_QA]->SetTitle("");
//      h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
//      h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
//      h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
//      h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
//      h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
//      h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
//
//      h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
//      h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
//      h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
//      h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
//      h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
//      h_mMass[KEY_QA]->GetYaxis()->CenterTitle();
//
//      h_mMass[KEY_QA]->SetMarkerStyle(24);
//      h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
//      h_mMass[KEY_QA]->SetMarkerSize(1.2);
//      h_mMass[KEY_QA]->Draw("pE");
//      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
//    }
//    if(i_theta == 10)
//    {
//      string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",ptQA,centQA,dcaQA,nsigQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
//
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Yields");
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();
//
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
//      h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
//      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
//      TF1 *f_bw = new TF1("f_bw",Poly3MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
//      f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][0]);
//      f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][1]);
//      f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][2]);
//      f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA][3]);
//      f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA][4]);
//      f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA][5]);
//      f_bw->SetLineColor(2);
//      f_bw->SetLineStyle(2);
//      f_bw->SetLineWidth(2);
//      f_bw->Draw("l same");
//
//      TF1 *f_bg = new TF1("f_bg",Poly3M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
//      f_bg->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][3]);
//      f_bg->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][4]);
//      f_bg->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][5]);
//      f_bg->SetLineColor(kBlack);
//      f_bg->SetLineStyle(2);
//      f_bg->SetLineWidth(2);
//      f_bg->Draw("l same");
//    }
//  }
//#endif

  /*
  // fit theta-differential bin to extract fit parameters for integration
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	  {
	    for(int i_theta = 0; i_theta < 10; i_theta++) // cos(theta*) loop
	    {
	      string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_Cos2PhiStarPhi_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	    }
	  }
	}
      }
    }
  }
  */

  // extract counts vs. pT with diffenretial integration ranges and methods
  vecFMap Par;
  TH1FMap h_mCounts;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;

  for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
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
                string KEY_rho = Form("rhoRaw_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                g_mRho[KEY_rho] = new TGraphAsymmErrors();
                for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; ++i_cent) // Centrality loop
                {
    	          string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  h_mCounts[KEY_counts] = new TH1F(KEY_counts.c_str(),KEY_counts.c_str(),10,-1.0,1.0);

                  string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  string KEY_InteTheta_Poly = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                  bool skip_y = false;
                   
                  g_mChiNDF[KEY_InteTheta_Poly] = new TGraphAsymmErrors();             
                  g_mPValue[KEY_InteTheta_Poly] = new TGraphAsymmErrors();             

                  for(int i_theta = 0; i_theta < 10; ++i_theta) // cos(theta*) loop
                  {
                    cout << "cent = " << i_cent << " i_dca = " << i_dca << " i_sig = " << i_sig << " i_norm = " << i_norm << " i_sigma = " << i_sigma << " i_method = " << i_method << " i_poly = " << i_poly << " i_pt = " << i_pt << " i_theta = " << i_theta << endl;
                    string KEY = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                    string KEY_Poly = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                    TF1 *f_bw;
                    if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
                    if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
                    if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
                    f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta_Poly][0]);
                    f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta_Poly][1]);
                    f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_Poly][2]/10.0);
                    f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_Poly][3]/10.0);
                    f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_Poly][4]/10.0);
                    if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_Poly][5]);
                    if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_Poly][6]);
                    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
                    TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"MQNRS");
                  
	            TF1 *f_bg;
                    if(i_poly == 0) f_bg = new TF1("f_bg",Poly, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
                    if(i_poly == 1) f_bg = new TF1("f_bg",Poly2,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
                    if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
                    TMatrixDSym covArr1(2);
                    TMatrixDSym covArr2(3);
                    TMatrixDSym covArr3(4);
                    double params1[2];
                    double params2[3];
                    double params3[4];
                    if(result.Get() != nullptr)// continue;
                    { 
                      f_bg->SetParameter(0,f_bw->GetParameter(3));
                      f_bg->SetParameter(1,f_bw->GetParameter(4));
                      if(i_poly >= 1) f_bg->SetParameter(2,f_bw->GetParameter(5));
                      if(i_poly >= 2) f_bg->SetParameter(3,f_bw->GetParameter(6));
                      f_bg->SetParError(0,f_bw->GetParError(3));
                      f_bg->SetParError(1,f_bw->GetParError(4));
                      if(i_poly >= 1) f_bg->SetParError(2,f_bw->GetParError(5));
                      if(i_poly >= 2) f_bg->SetParError(3,f_bw->GetParError(6));
                  
                      if(i_poly == 0) 
                      {
                        params1[0] = result->GetParams()[3];
                        params1[1] = result->GetParams()[4];
                        covArr1(0,0) = result->GetCovarianceMatrix()(3,3);
                        covArr1(0,1) = result->GetCovarianceMatrix()(3,4);
                        covArr1(1,0) = result->GetCovarianceMatrix()(4,3);
                        covArr1(1,1) = result->GetCovarianceMatrix()(4,4);
                      }
                      if(i_poly == 1) 
                      {
                        params2[0] = result->GetParams()[3];
                        params2[1] = result->GetParams()[4];
                        params2[2] = result->GetParams()[5];
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
                        params3[0] = result->GetParams()[3];
                        params3[1] = result->GetParams()[4];
                        params3[2] = result->GetParams()[5];
                        params3[3] = result->GetParams()[6];
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
                    }
                    else
  	            { 
                      skip_y = true;  

	              f_bg->SetParameter(0,0.0);
	              f_bg->SetParameter(1,0.0);
	              if(i_poly >= 1) f_bg->SetParameter(2,0.0);
	              if(i_poly >= 2) f_bg->SetParameter(3,0.0);
	              f_bg->SetParError(0,1.0);
	              f_bg->SetParError(1,1.0);
	              if(i_poly >= 1) f_bg->SetParError(2,1.0);
	              if(i_poly >= 2) f_bg->SetParError(3,1.0);
	         
                      if(i_poly == 0)
                      { 
                        params1[0] = 0.0;//result->GetParams()[3];
                        params1[1] = 0.0;//result->GetParams()[4];
                        covArr1(0,0) = 1.0;//result->GetCovarianceMatrix()(3,3);
                        covArr1(0,1) = 0.0;//result->GetCovarianceMatrix()(3,4);
                        covArr1(1,0) = 0.0;//result->GetCovarianceMatrix()(4,3);
                        covArr1(1,1) = 1.0;//result->GetCovarianceMatrix()(4,4);
                      }
                      if(i_poly == 1)
                      { 
                        params2[0] = 0.0;//result->GetParams()[3];
                        params2[1] = 0.0;//result->GetParams()[4];
                        params2[2] = 0.0;//result->GetParams()[5];
                        covArr2(0,0) = 1.0;//result->GetCovarianceMatrix()(3,3);
                        covArr2(0,1) = 0.0;//result->GetCovarianceMatrix()(3,4);
                        covArr2(0,2) = 0.0;//result->GetCovarianceMatrix()(3,5);
                        covArr2(1,0) = 0.0;//result->GetCovarianceMatrix()(4,3);
                        covArr2(1,1) = 1.0;//result->GetCovarianceMatrix()(4,4);
                        covArr2(1,2) = 0.0;//result->GetCovarianceMatrix()(4,5);
                        covArr2(2,0) = 0.0;//result->GetCovarianceMatrix()(5,3);
                        covArr2(2,1) = 0.0;//result->GetCovarianceMatrix()(5,4);
                        covArr2(2,2) = 1.0;//result->GetCovarianceMatrix()(5,5);
                      }
                      if(i_poly == 2)
                      { 
                        params3[0] = 0.0;//result->GetParams()[3];
                        params3[1] = 0.0;//result->GetParams()[4];
                        params3[2] = 0.0;//result->GetParams()[5];
                        params3[3] = 0.0;//result->GetParams()[5];
                        covArr3(0,0) = 1.0;//result->GetCovarianceMatrix()(3,3);
                        covArr3(0,1) = 0.0;//result->GetCovarianceMatrix()(3,4);
                        covArr3(0,2) = 0.0;//result->GetCovarianceMatrix()(3,5);
                        covArr3(0,3) = 0.0;//result->GetCovarianceMatrix()(3,6);
                        covArr3(1,0) = 0.0;//result->GetCovarianceMatrix()(4,3);
                        covArr3(1,1) = 1.0;//result->GetCovarianceMatrix()(4,4);
                        covArr3(1,2) = 0.0;//result->GetCovarianceMatrix()(4,5);
                        covArr3(1,3) = 0.0;//result->GetCovarianceMatrix()(4,6);
                        covArr3(2,0) = 0.0;//result->GetCovarianceMatrix()(5,3);
                        covArr3(2,1) = 0.0;//result->GetCovarianceMatrix()(5,4);
                        covArr3(2,2) = 1.0;//result->GetCovarianceMatrix()(5,5);
                        covArr3(2,3) = 0.0;//result->GetCovarianceMatrix()(5,6);
                        covArr3(3,0) = 0.0;//result->GetCovarianceMatrix()(6,3);
                        covArr3(3,1) = 0.0;//result->GetCovarianceMatrix()(6,4);
                        covArr3(3,2) = 0.0;//result->GetCovarianceMatrix()(6,5);
                        covArr3(3,3) = 1.0;//result->GetCovarianceMatrix()(6,6);
                      }
 
                      f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta_Poly][0]);
                      f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta_Poly][1]);
                      f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_Poly][2]/10.0);
                      f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_Poly][3]/10.0);
                      f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_Poly][4]/10.0);
                      if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_Poly][5]);
                      if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_Poly][6]);
                      f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	            }
                    
                                     

                    float bin_width = h_mMass[KEY]->GetBinWidth(1);
                    float Inte_start = Par_InteTheta[KEY_InteTheta_Poly][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]-0.5*bin_width;
                    float Inte_stop  = Par_InteTheta[KEY_InteTheta_Poly][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]+0.5*bin_width;
                    float counts_bg;
                    if(result.Get() != nullptr) counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
                    if(result.Get() == nullptr) counts_bg = 0.0;// continue;
                    float errors_bg;
                    if(result.Get() != nullptr) 
                    {
                      if(i_poly == 0) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params1,covArr1.GetMatrixArray())/bin_width;
                      if(i_poly == 1) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params2,covArr2.GetMatrixArray())/bin_width;
                      if(i_poly == 2) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params3,covArr3.GetMatrixArray())/bin_width;
                    }
                    if(result.Get() == nullptr) errors_bg = 0.0;// continue;
                  
	            float cos_mean = (float(i_theta)-4.5)/5.0;
                    float chi2NDF = float(f_bw->GetChisquare())/float(f_bw->GetNDF()); 
                    float pvalue = TMath::Prob(f_bw->GetChisquare(),f_bw->GetNDF());
                    g_mChiNDF[KEY_InteTheta_Poly]->SetPoint(i_theta,cos_mean,chi2NDF);
                    g_mPValue[KEY_InteTheta_Poly]->SetPoint(i_theta,cos_mean,pvalue);
                    float bin_center = (float(i_theta)-4.5)/5.0;
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
                      h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts-counts_bg);
                      h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors+errors_bg*errors_bg));
                    }
                    if(i_method == 1)
                    {
                      float counts_bw;
                      if(result.Get() != nullptr) counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
                      if(result.Get() == nullptr) counts_bw = 0.0;// continue;
                      float errors_bw; 
                      if(result.Get() != nullptr) errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
                      if(result.Get() == nullptr) errors_bw = 0.0;// continue;
                      h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts_bw-counts_bg);
                      h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));
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
                    //TF1 *f_bg = new TF1("f_bg",Poly3M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
                    //f_bg->SetParameter(0,f_bw->GetParameter(3));
                    //f_bg->SetParameter(1,f_bw->GetParameter(4));
                    //f_bg->SetParameter(2,f_bw->GetParameter(5));
                    //f_bg->SetParError(0,f_bw->GetParError(3));
                    //f_bg->SetParError(1,f_bw->GetParError(4));
                    //f_bg->SetParError(2,f_bw->GetParError(5));
                  
                    //double params[3] = {result->GetParams()[3],result->GetParams()[4],result->GetParams()[5]};
                    //TMatrixDSym covArr(3);
                    //covArr(0,0) = result->GetCovarianceMatrix()(3,3);
                    //covArr(0,1) = result->GetCovarianceMatrix()(3,4);
                    //covArr(0,2) = result->GetCovarianceMatrix()(3,5);
                    //covArr(1,0) = result->GetCovarianceMatrix()(4,3);
                    //covArr(1,1) = result->GetCovarianceMatrix()(4,4);
                    //covArr(1,2) = result->GetCovarianceMatrix()(4,5);
                    //covArr(2,0) = result->GetCovarianceMatrix()(5,3);
                    //covArr(2,1) = result->GetCovarianceMatrix()(5,4);
                    //covArr(2,2) = result->GetCovarianceMatrix()(5,5);
                  
                    //float bin_width = h_mMass[KEY]->GetBinWidth(1);
                    //float Inte_start = Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
                    //float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
                    //float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
                    //float errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params,covArr.GetMatrixArray())/bin_width;
                  
                    //float bin_center = 1/14.0+i_theta/10.0;
                    //if(i_method == 0)
                    //{
                    //  int bin_start = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
                    //  int bin_stop  = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
                    //  float counts = 0.0;
                    //  float errors = 0.0;
                    //  for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
                    //  {
                    //    counts += h_mMass[KEY]->GetBinContent(i_bin);
                    //    errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
                    //  }
                    //  h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts-counts_bg);
                    //  h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors+errors_bg*errors_bg));
                    //}
                    //if(i_method == 1)
                    //{
                    //  float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
                    //  float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
                    //  h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts_bw-counts_bg);
                    //  h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));
                    //}
                    //Par[KEY].clear();
                    //Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(0)));
                    //Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(1)));
                    //Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(2)));
                    //Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(3)));
                    //Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(4)));
                    //Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(5)));
                 
                  //float cent_mean = (vmsa::centVal[i_cent]+vmsa::centVal[i_cent+1])/2.0;
                  float cent_mean = (vmsa::centVal[i_cent]+vmsa::centVal[i_cent+1])/2.0; 
                  TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
                  f_rho->SetParameter(0,0.33);
                  f_rho->SetParameter(1,h_mCounts[KEY_counts]->GetMaximum());
                  h_mCounts[KEY_counts]->Fit(f_rho,"NMRI");
                  Par_rhoFit[KEY_counts].clear();
                  //Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
                  //Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
                  //g_mRho[KEY_rho]->SetPoint(i_pt,cent_mean,f_rho->GetParameter(0));
                  //g_mRho[KEY_rho]->SetPointError(i_pt,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
                  if(!skip_y)
                  {
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
                    g_mRho[KEY_rho]->SetPoint(i_cent,cent_mean,f_rho->GetParameter(0));
                    g_mRho[KEY_rho]->SetPointError(i_cent,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
                  }
                  else
                  {
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(-999.0));
                    Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
                    g_mRho[KEY_rho]->SetPoint(i_cent,cent_mean,-999.0);
                    g_mRho[KEY_rho]->SetPointError(i_cent,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
 
  cout << "Created rho00 plots" << endl;


#if _PlotQA_
  for(int i_poly = 0; i_poly < 3; i_poly++)
  {
    string outputname = Form("./figures/%s/%s/centralitystudy/Cos2PhiStarPhi_allThetaYields%s_%s_Order%d_%s_CentDependence_Poly%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str(),i_poly+1);
    string output_start = Form("%s[",outputname.c_str());
    TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
    c_diff->Print(output_start.c_str());
    for(int i_cent = vmsa::Cent_start; i_cent < 9; ++i_cent) // Centrality loop
    {  
      //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
      for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
      {
        c_diff->Clear();
        c_diff->Divide(4,3);
        for(int i_theta = 0; i_theta < 12; ++i_theta)
        {
          c_diff->cd(i_theta+1);
          c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
          c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
          c_diff->cd(i_theta+1)->SetTicks(1,1);
          c_diff->cd(i_theta+1)->SetGrid(0,0);
          if(i_theta < 10)
          {
            string KEY_QA = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
            h_mMass[KEY_QA]->SetTitle("");
            h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
            h_mMass[KEY_QA]->SetTitle(Form("%.2f<p_{T}<%.2f, %d/7<|cos(#theta*)|<%d/7, cent%d",vmsa::pt_low_cent[energy][i_pt],vmsa::pt_up_cent[energy][i_pt],i_theta,i_theta+1,i_cent));
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

            //string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
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
          if(i_theta == 10)
          {
            string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
            string KEY_InteTheta_QA_Poly = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,i_poly+1);
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle(Form("%.2f<p_{T}<%.2f, 0/7<|cos(#theta*)|<7/7, cent%d",vmsa::pt_low_cent[energy][i_pt],vmsa::pt_up_cent[energy][i_pt],i_cent));
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
            TF1 *f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
            f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA_Poly][0]);
            f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA_Poly][1]);
            f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA_Poly][2]);
            f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA_Poly][3]);
            f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA_Poly][4]);
            f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA_Poly][5]);
            f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_QA_Poly][6]);
            f_bw->SetLineColor(2);
            f_bw->SetLineStyle(2);
            f_bw->SetLineWidth(2);
            f_bw->Draw("l same");
          }
        }

        for(int i_theta = 0; i_theta < 10; ++i_theta)
        {
          c_diff->cd(i_theta+1);
          string KEY_QA = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
          string KEY_QA_Poly = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,i_theta,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,i_poly+1);
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

        c_diff->cd(11);
        for(int i_sigma = vmsa::Sig_start; i_sigma < /*vmsa::Sig_stop*/1; ++i_sigma)
        {
          for(int i_method = /*vmsa::Method_start*/1; i_method < vmsa::Method_stop; ++i_method)
          {
            string KEY_counts_QA = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            if(i_sigma == vmsa::Sig_start && i_method == vmsa::Method_start)
            {
              h_mCounts[KEY_counts_QA];
              h_mCounts[KEY_counts_QA]->SetTitle("");
              h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
              h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
              h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();

              h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
              h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("Counts");
              h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();

              h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
              h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
              h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
              h_mCounts[KEY_counts_QA]->Draw("pE");
            }
            else
            {
              h_mCounts[KEY_counts_QA]->SetMarkerStyle(24);
              h_mCounts[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
              h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
              h_mCounts[KEY_counts_QA]->Draw("pE same");
            }
            TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
            f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
            f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
            f_rho->SetLineColor(i_sigma+10*i_method+1);
            f_rho->SetLineWidth(2);
            f_rho->SetLineStyle(2);
            f_rho->Draw("l same");
          }
        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());
      }
    }
    string output_stop = Form("%s]",outputname.c_str());
    c_diff->Print(output_stop.c_str()); // close pdf file
   // if(!random3D) c_diff->SaveAs("../figures/c_diff_2.pdf");
   // if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2.pdf");
  }
#endif

    
//#if _PlotQA_
//  for(int i_theta = 0; i_theta < 10; ++i_theta)
//  {
//    c_diff->cd(i_theta+1);
//    string KEY_QA = Form("pt_%d_Centrality_%d_Cos2PhiStarPhi_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d",ptQA,centQA,i_theta,dcaQA,nsigQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
//  //  cout << KEY_QA << endl;
//    TF1 *f_bw = new TF1("f_bw",Poly3MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
//    f_bw->SetParameter(0,Par[KEY_QA][0]);
//    f_bw->SetParameter(1,Par[KEY_QA][1]);
//    f_bw->SetParameter(2,Par[KEY_QA][2]);
//    f_bw->SetParameter(3,Par[KEY_QA][3]);
//    f_bw->SetParameter(4,Par[KEY_QA][4]);
//    f_bw->SetParameter(5,Par[KEY_QA][5]);
//    f_bw->SetLineColor(4);
//    f_bw->SetLineStyle(2);
//    f_bw->SetLineWidth(2);
//    f_bw->Draw("l same");
//
//    TF1 *f_bg = new TF1("f_bg",Poly3M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
//    f_bg->SetParameter(0,Par[KEY_QA][3]);
//    f_bg->SetParameter(1,Par[KEY_QA][4]);
//    f_bg->SetParameter(2,Par[KEY_QA][5]);
//    f_bg->SetLineColor(kRed);
//    f_bg->SetLineStyle(2);
//    f_bg->SetLineWidth(2);
//    f_bg->Draw("l same");
//  }
//
//  c_diff->cd(9);
//  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
//  {
//    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
//    {
//      string KEY_counts_QA = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ptQA,centQA,dcaQA,nsigQA,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,i_sigma,vmsa::mInteMethod[i_method].c_str());
//      //cout << KEY_counts_QA << endl;
//      if(i_sigma == vmsa::Sig_start && i_method == vmsa::Method_start)
//      {
//	h_mCounts[KEY_counts_QA];
//	h_mCounts[KEY_counts_QA]->SetTitle("");
//	h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
//	h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
//	h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta*)");
//	h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
//	h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
//	h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();
//
//	h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
//	h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
//	h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("Counts");
//	h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
//	h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
//	h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();
//
//	h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
//	h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
//	h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
//	h_mCounts[KEY_counts_QA]->Draw("pE");
//      }
//      else
//      {
//	h_mCounts[KEY_counts_QA]->SetMarkerStyle(24);
//	h_mCounts[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
//	h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
//	h_mCounts[KEY_counts_QA]->Draw("pE same");
//      }
//      TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
//      f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
//      f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
//      f_rho->SetLineColor(i_sigma+10*i_method+1);
//      f_rho->SetLineWidth(2);
//      f_rho->SetLineStyle(2);
//      f_rho->Draw("l same");
//    }
//  }
//  if(!random3D) c_diff->SaveAs("../figures/c_diff_2_poly1_phi.pdf");
//  if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2.pdf");
//#endif

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,10.05);
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

  h_frame->GetYaxis()->SetRangeUser(0.25,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();

  cout << "Created h_frame" << endl;
  h_frame->DrawCopy("pE");
  PlotLine(0.0,80.0,1.0/3.0,1.0/3.0,1,2,2);

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
            string KEY_rho = Form("rhoRaw_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly3",ptQA,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
  //          cout << KEY_rho << endl;
            Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,0,1.1);
          }
        }
      }
    }
  }
  if(!random3D) c_rho->SaveAs("../figures/c_rho.pdf");
  if(random3D)  c_rho->SaveAs("../figures/3DRandom/c_rho.pdf");

  TCanvas *c_chi = new TCanvas("c_chi","c_chi",10,10,1200,1200);
  c_chi->Divide(3,3);
  for(int i = 0; i < 9; i++)
  {
    c_chi->cd(i+1);
    c_chi->cd(i+1)->SetLeftMargin(0.15);
    c_chi->cd(i+1)->SetBottomMargin(0.15);
    c_chi->cd(i+1)->SetTicks(1,1);
    c_chi->cd(i+1)->SetGrid(0,0);
  }  

  int polycolor[3] = {kGray+2,kBlue,kOrange+7};
  int polystyle[3] = {24,25,26};

  string outputname = Form("figures/%s/%s/centralitystudy/Cos2PhiStarPhi_Chi2NDF_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  string output_start = Form("%s[",outputname.c_str());
  
  c_chi->Print(output_start.c_str());
  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      if( i_dca != 0 && i_sig != 0 ) continue;
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
        for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
        {
          string Title;           
          TLegend *leg1 = new TLegend(0.75,0.75,0.9,0.9);
          for(int i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
          {
            c_chi->cd(i_cent+1);
            double min[3], max[3];
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_pT_Poly = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
              min[i_poly] = TMath::MinElement(7,g_mChiNDF[KEY_pT_Poly]->GetY());
              max[i_poly] = TMath::MaxElement(7,g_mChiNDF[KEY_pT_Poly]->GetY());
              cout << "Min = " << min[i_poly] << ", Max = " << max[i_poly] << endl;
              g_mChiNDF[KEY_pT_Poly]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
              g_mChiNDF[KEY_pT_Poly]->GetYaxis()->SetTitle("#chi^{2}/NDF");
              g_mChiNDF[KEY_pT_Poly]->SetMarkerColor(polycolor[i_poly]);
              g_mChiNDF[KEY_pT_Poly]->SetMarkerStyle(polystyle[i_poly]);
              if(i_cent == 0) leg1->AddEntry(g_mChiNDF[KEY_pT_Poly],Form("Poly %d",i_poly+1),"p");
            }
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_pT_Poly = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
              if(i_poly == 0)
              {
                g_mChiNDF[KEY_pT_Poly]->SetMinimum(TMath::MinElement(3,min)*0.8);
                g_mChiNDF[KEY_pT_Poly]->SetMaximum(TMath::MaxElement(3,max)*1.4);
                g_mChiNDF[KEY_pT_Poly]->SetTitle(KEY_pT_Poly.c_str());
                g_mChiNDF[KEY_pT_Poly]->Draw("AP");
              }
              else 
              {
                g_mChiNDF[KEY_pT_Poly]->Draw("P same");
              }
            }
            leg1->Draw("same");
          }

          //TLegend *leg2 = new TLegend(0.2,0.2,0.4,0.4);
          //c_chi->cd(5);
          //Title = Form("Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);

          //for(int i_poly = 0; i_poly < 3; i_poly++)
          //{
          //  string KEY_Poly = Form("Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //  g_mChiNDF[KEY_Poly]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          //  g_mChiNDF[KEY_Poly]->GetYaxis()->SetTitle("#chi^{2}/NDF");
          //  g_mChiNDF[KEY_Poly]->SetMarkerColor(polycolor[i_poly]);
          //  g_mChiNDF[KEY_Poly]->SetMarkerStyle(polystyle[i_poly]);
          //  if(i_poly == 0)
          //  {
          //    g_mChiNDF[KEY_Poly]->SetTitle(Title.c_str());
          //    g_mChiNDF[KEY_Poly]->Draw("AP");
          //  }
          //  else 
          //  {
          //    g_mChiNDF[KEY_Poly]->Draw("P same");
          //  }
          //  leg2->AddEntry(g_mChiNDF[KEY_Poly],Form("Poly %d",i_poly),"p");
          //}
          //leg2->Draw("same");
          c_chi->Update();
          c_chi->Print(outputname.c_str());
          
        }
      }
    }  
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_chi->Print(output_stop.c_str()); // close pdf file

  outputname = Form("figures/%s/%s/centralitystudy/Cos2PhiStarPhi_PValue_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  output_start = Form("%s[",outputname.c_str());
  
  c_chi->Print(output_start.c_str());

  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      if( i_dca != 0 && i_sig != 0 ) continue;
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
        for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
        {
          string Title;           
          TLegend *leg1 = new TLegend(0.2,0.2,0.4,0.4);
          for(int i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
          {
            c_chi->cd(i_cent+1);
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_pT_Poly = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
              g_mPValue[KEY_pT_Poly]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
              g_mPValue[KEY_pT_Poly]->GetYaxis()->SetTitle("p-value");
              g_mPValue[KEY_pT_Poly]->SetMarkerColor(polycolor[i_poly]);
              g_mPValue[KEY_pT_Poly]->SetMarkerStyle(polystyle[i_poly]);
              if(i_poly == 0)
              {
                g_mPValue[KEY_pT_Poly]->SetTitle(KEY_pT_Poly.c_str());
                g_mPValue[KEY_pT_Poly]->Draw("AP");
              }
              else 
              {
                g_mPValue[KEY_pT_Poly]->Draw("P same");
              }
              if(i_cent == 0) leg1->AddEntry(g_mPValue[KEY_pT_Poly],Form("Poly %d",i_poly+1),"p");
            }
            leg1->Draw("same");
          }

          //TLegend *leg2 = new TLegend(0.2,0.2,0.4,0.4);
          //c_chi->cd(5);
          //Title = Form("Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);

          //for(int i_poly = 0; i_poly < 3; i_poly++)
          //{
          //  string KEY_Poly = Form("Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //  g_mPValue[KEY_Poly]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          //  g_mPValue[KEY_Poly]->GetYaxis()->SetTitle("p-value");
          //  g_mPValue[KEY_Poly]->SetMarkerColor(polycolor[i_poly]);
          //  g_mPValue[KEY_Poly]->SetMarkerStyle(polystyle[i_poly]);
          //  if(i_poly == 0)
          //  {
          //    g_mPValue[KEY_Poly]->SetTitle(Title.c_str());
          //    g_mPValue[KEY_Poly]->Draw("AP");
          //  }
          //  else 
          //  {
          //    g_mPValue[KEY_Poly]->Draw("P same");
          //  }
          //  leg2->AddEntry(g_mPValue[KEY_Poly],Form("Poly %d",i_poly),"p");
          //}
          //leg2->Draw("same");
          c_chi->Update();
          c_chi->Print(outputname.c_str());
        }
      }
    }
  }
  output_stop = Form("%s]",outputname.c_str());
  c_chi->Print(output_stop.c_str()); // close pdf file

  string outputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_RawRhoCentSys_%s_Poly.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_RawRhoCentSys_%s_Poly_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/RawRhoCentSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
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
                string KEY_rho = Form("rhoRaw_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
                g_mRho[KEY_rho]->Write();
                //cout << KEY_rho << endl;
                for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; ++i_cent) // Centrality loop
                {
                  string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  h_mCounts[KEY_counts]->Write();
                  //cout << KEY_counts << endl;
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
