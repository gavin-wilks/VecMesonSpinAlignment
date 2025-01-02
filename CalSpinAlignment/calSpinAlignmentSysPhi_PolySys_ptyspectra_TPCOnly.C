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


void calSpinAlignmentSysPhi_PolySys_ptyspectra_TPCOnly(int energy = 4, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1")
{

  gStyle->SetNumberContours(250);
  
  std::string EP[2] = {"","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/TPCOnly_ptyspectra_InvMassSubBg_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/InvMassSubBg_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
 
  TGraMap g_mChiNDF;
  TGraMap g_mPValue;
  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  for(int i_pt = 0; i_pt < vmsa::rebinpttotal; i_pt++) // pt loop
  {
    for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    for(int i_eta = 0; i_eta < vmsa::rebinytotal; i_eta++) // cos(theta*) loop
	    {
	      string KEY = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              cout << KEY << endl;
	      h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str());
            }
	  }
	}
      }
    }
  }

  const int nybins = 24;
  double ybins[nybins+1] = {0.0};
 
  bool skip[10][25][24] = {false};

  for(int iy = 0; iy < nybins+1; iy++)
  {
    ybins[iy] = double(iy-12.)/8.;
    cout << "iy = " << iy << ", y = " << ybins[iy] << endl;
  }

  // extract counts vs. pT with diffenretial integration ranges and methods
  float InvMassLow[vmsa::rebinpttotal][vmsa::rebinytotal] = {0.0};
  float InvMassHigh[vmsa::rebinpttotal][vmsa::rebinytotal] = {0.0};

  float fitparams[vmsa::rebinpttotal][vmsa::rebinytotal][5] = {0.0};

  float fpoly[vmsa::rebinpttotal][vmsa::rebinytotal] = {0.0};

  vecFMap Par;
  TH1FMap h_mCounts;
  TH2FMap h_mSpectra;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;
  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
              for(int i_poly = 0; i_poly < 3; i_poly++)
              {
	        string KEY_spectra = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                cout << KEY_spectra << endl;
	        //h_mSpectra[KEY_spectra] = new TH2F(KEY_spectra.c_str(),KEY_spectra.c_str(),vmsa::pt_total,vmsa::pt_bin,nybins,ybins);
	        h_mSpectra[KEY_spectra] = new TH2F(KEY_spectra.c_str(),KEY_spectra.c_str(),vmsa::rebinpttotal,vmsa::rebinptval,vmsa::rebinytotal,vmsa::rebinyval);
	        h_mSpectra[KEY_spectra]->SetStats(0);
                cout << "Created 2D plot" << endl;
	        for(int i_pt = 0/*vmsa::pt_start*/; i_pt < vmsa::rebinpttotal/*< vmsa::pt_stop*/; ++i_pt) // pt loop
	        {
                  cout << "i_pt = " << i_pt << endl;
	          for(int i_eta = 0; i_eta < vmsa::rebinytotal; ++i_eta) // cos(theta*) loop
	          {
                    cout << "i_eta = " << i_eta << endl;
                    double pt_mean = (vmsa::rebinptval[i_pt]+vmsa::rebinptval[i_pt+1])/2.0;
                    double y_mean =  (vmsa::rebinyval[i_eta]+vmsa::rebinyval[i_eta+1])/2.0;
                     
	            string KEY = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	            string KEY_Poly = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
                    cout << KEY << endl;
                    cout << KEY_Poly << endl;
	            TF1 *f_bw;
                    if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
                    if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
                    if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
	            f_bw->SetParameter(0,vmsa::InvMass[pid]);
	            f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.005,vmsa::InvMass[pid]+0.005);
	            //f_bw->SetParameter(1,vmsa::Width[pid]);
	            f_bw->SetParameter(1,0.006);
	            f_bw->SetParameter(2,1.0);
	            float norm = h_mMass[KEY]->GetMaximum();///f_bw->GetMaximum();
	            f_bw->SetParameter(2,norm);
	            f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
                    //if(h_mMass[KEY]->GetEntries() <= 2) continue;
                    //cout << h_mMass[KEY]->GetEntries() << endl;
	            TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"MQNRS");
                    if(result.Get() == nullptr) 
                    {
                      skip[i_cent][i_pt][i_eta] = true;
                      cout << "Result is NULL SKIP!!!!" << endl;
                      continue;
                    }
                    //cout << "Just fit the result" << endl;
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
	            float Inte_start = f_bw->GetParameter(0)-vmsa::nSigVecSys[i_sigma]*f_bw->GetParameter(1)-0.5*bin_width;
                    float Inte_stop  = f_bw->GetParameter(0)+vmsa::nSigVecSys[i_sigma]*f_bw->GetParameter(1)+0.5*bin_width;
                    
                    InvMassLow[i_pt][i_eta]  = Inte_start;
                    InvMassHigh[i_pt][i_eta] = Inte_stop;
 
                    float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
                    float errors_bg; 
                    if(i_poly == 0) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params1,covArr1.GetMatrixArray())/bin_width;
                    if(i_poly == 1) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params2,covArr2.GetMatrixArray())/bin_width;
                    if(i_poly == 2) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params3,covArr3.GetMatrixArray())/bin_width;

	            if(i_method == 0)
	            {
	              int bin_start = h_mMass[KEY]->FindBin(f_bw->GetParameter(0)-vmsa::nSigVecSys[i_sigma]*f_bw->GetParameter(1));
	              int bin_stop  = h_mMass[KEY]->FindBin(f_bw->GetParameter(0)+vmsa::nSigVecSys[i_sigma]*f_bw->GetParameter(1));
	              float counts = 0.0;
	              float errors = 0.0;
	              for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	              {
	                counts += h_mMass[KEY]->GetBinContent(i_bin);
	                errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	              }
                      float final_counts = counts-counts_bg;
                      if(final_counts < 0.0) final_counts = 0.0;
	              h_mSpectra[KEY_spectra]->SetBinContent(h_mSpectra[KEY_spectra]->FindBin(pt_mean,y_mean),final_counts);
	              h_mSpectra[KEY_spectra]->SetBinError(h_mSpectra[KEY_spectra]->FindBin(pt_mean,y_mean),TMath::Sqrt(errors+errors_bg*errors_bg));
	            }
	            if(i_method == 1)
	            {
	              //float bin_width = h_mMass[KEY]->GetBinWidth(1);
	              //float Inte_start = Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
	              //float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
	              float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
	              float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
                      float final_counts_bw = counts_bw-counts_bg;
                      if(final_counts_bw < 0.0) final_counts_bw = 0.0; 
	              h_mSpectra[KEY_spectra]->SetBinContent(h_mSpectra[KEY_spectra]->FindBin(pt_mean,y_mean),final_counts_bw);
	              h_mSpectra[KEY_spectra]->SetBinError(h_mSpectra[KEY_spectra]->FindBin(pt_mean,y_mean),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));

                      fpoly[i_pt][i_eta] = counts_bg/counts_bw;
	            }
	            Par[KEY_Poly].clear();
	            Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(0)));
                    //cout << "Push back 0: " << Par[KEY_Poly][0] << endl;
	            Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(1)));
                    //cout << "Push back 1: " << Par[KEY_Poly][1] << endl;
	            Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(2)));
                    //cout << "Push back 2: " << Par[KEY_Poly][2] << endl;
	            Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(3)));
                    //cout << "Push back 3: " << Par[KEY_Poly][3] << endl;
	            Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(4)));
                    
                    //cout << "Push back 4: " << Par[KEY_Poly][4] << endl;
                    for(int i = 0; i < 5; i++) fitparams[i_pt][i_eta][i] = f_bw->GetParameter(i);
                     
	            if(i_poly >= 1) Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(5)));
	            if(i_poly >= 2) Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(6)));
	          }
	        }
              }
	    }
	  }
	}
      }
    }
  }

  cout << endl << "const float InvMassLow[rebinpttotal][rebinytotal] = {";
  for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
  {
    cout << "{";
    for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
    {
      if(etabin < vmsa::rebinytotal-1) cout << InvMassLow[ptbin][etabin] << ",";
      else                             cout << InvMassLow[ptbin][etabin];
    }
    if(ptbin < vmsa::rebinpttotal -1 ) cout << "}," << endl;
    else                                cout << "}";
  }
  cout << "};" << endl << endl;

  cout << endl << "const float InvMassHigh[rebinpttotal][rebinytotal] = {";
  for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
  {
    cout << "{";
    for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
    {
      if(etabin < vmsa::rebinytotal-1) cout << InvMassHigh[ptbin][etabin] << ",";
      else                             cout << InvMassHigh[ptbin][etabin];
    }
    if(ptbin < vmsa::rebinpttotal -1 ) cout << "}," << endl;
    else                                cout << "}";
  }
  cout << "};" << endl << endl;

  cout << endl << "const float fitparams[rebinpttotal][rebinytotal][5] = {";
  for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
  {
    cout << "{";
    for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
    {
      cout << "{";
      for(int i  = 0; i < 5; i++)
      {
        if(i < 4) cout << fitparams[ptbin][etabin][i] << ",";
        else      cout << fitparams[ptbin][etabin][i];
      }
      if(etabin < vmsa::rebinytotal-1) cout << "}," << endl;
      else                             cout << "}";
    }
    if(ptbin < vmsa::rebinpttotal -1 ) cout << "}," << endl;
    else                               cout << "}";
  }
  cout << "};" << endl << endl;



  cout << endl << "const float fpoly[rebinpttotal][rebinytotal] = {";
  for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
  {
    cout << "{";
    for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
    {
      if(etabin < vmsa::rebinytotal-1) cout << fpoly[ptbin][etabin] << ",";
      else                             cout << fpoly[ptbin][etabin];
    }
    if(ptbin < vmsa::rebinpttotal -1 ) cout << "}," << endl;
    else                                cout << "}";
  }
  cout << "};" << endl << endl;

#if _PlotQA_
  for(int i_poly = 0; i_poly < 3; i_poly++)
  {
    for(int i_norm = 0; i_norm < 3; i_norm++)
    {
      for(int i_cent = 9; i_cent < 10; ++i_cent) // Centrality loop
      {
        string outputname = Form("./figures/%s/%s/pTstudy/TPCOnly_ptyspectra_allThetaYields%s_%s_Order%d_cent%d_%s_pTdependence_Norm%d_Poly%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_cent,etamode.c_str(),i_norm,i_poly+1);
        string output_start = Form("%s[",outputname.c_str());
        
        //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
        TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,2400,1600);
        c_diff->Print(output_start.c_str());
        //for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	for(int i_pt = 0/*vmsa::pt_start*/; i_pt < vmsa::rebinpttotal/*< vmsa::pt_stop*/; ++i_pt) // pt loop
        {
          c_diff->Clear();
          c_diff->Divide(6,4);
	  for(int i_eta = 0; i_eta < vmsa::rebinytotal; i_eta++) // cos(theta*) loop
	  {
            c_diff->cd(i_eta+1);
            c_diff->cd(i_eta+1)->SetLeftMargin(0.15);
            c_diff->cd(i_eta+1)->SetBottomMargin(0.15);
            c_diff->cd(i_eta+1)->SetTicks(1,1);
            c_diff->cd(i_eta+1)->SetGrid(0,0);
            if(skip[i_cent][i_pt][i_eta]) continue;
            string KEY_QA = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),i_norm);
            //cout << KEY_QA << endl;
            h_mMass[KEY_QA]->SetTitle("");
            h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
            h_mMass[KEY_QA]->SetTitle(Form("%.2f<p_{T}<%.2f, %.2f<y<%.2f",vmsa::rebinptval[i_pt],vmsa::rebinptval[i_pt+1],vmsa::rebinyval[i_eta],vmsa::rebinyval[i_eta+1]));
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
  
            string KEY_QA_Poly = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Poly%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
            //cout << KEY_QA_Poly << endl;
            TF1 *f_bw;
            if(i_poly == 0) f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
            if(i_poly == 1) f_bw = new TF1("f_bw",Poly2BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
            if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
            //cout << "About to print parameters" << endl;
            f_bw->SetParameter(0,Par[KEY_QA_Poly][0]);
            //cout << "Parameter 0 = " << Par[KEY_QA_Poly][0] << endl;
            f_bw->SetParameter(1,Par[KEY_QA_Poly][1]);
            //cout << "Parameter 1 = " << Par[KEY_QA_Poly][1] << endl;
            f_bw->SetParameter(2,Par[KEY_QA_Poly][2]);
            //cout << "Parameter 2 = " << Par[KEY_QA_Poly][2] << endl;
            f_bw->SetParameter(3,Par[KEY_QA_Poly][3]);
            //cout << "Parameter 3 = " << Par[KEY_QA_Poly][3] << endl;
            f_bw->SetParameter(4,Par[KEY_QA_Poly][4]);
            //cout << "Parameter 4 = " << Par[KEY_QA_Poly][4] << endl;
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
        string output_stop = Form("%s]",outputname.c_str());
        c_diff->Print(output_stop.c_str()); // close pdf file
      }
    }
     // if(!random3D) c_diff->SaveAs("../figures/c_diff_2.pdf");
     // if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2.pdf");
  }

  for(int i_poly = 0; i_poly < 3/*3*/; i_poly++)
  {
    for(int i_norm = 0; i_norm < 3; i_norm++)
    {
      for(int i_cent = 9; i_cent < 10; ++i_cent) // Centrality loop
      {
        string outputname = Form("./figures/%s/%s/pTstudy/TPCOnly_ptyspectra_2Dyield%s_%s_Order%d_cent%d_%s_pTdependence_Norm%d_Poly%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_cent,etamode.c_str(),i_norm,i_poly+1);
        
        //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
        TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,400,400);
        c_diff->SetLeftMargin(0.15);
        c_diff->SetBottomMargin(0.15);
        c_diff->SetRightMargin(0.15);
        c_diff->SetTicks(1,1);
        c_diff->SetGrid(0,0);
        c_diff->SetLogz();
        for(int i_sigma = 0/*vmsa::Sig_start*/; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
        {
          for(int i_method = 1/*vmsa::Method_start*/; i_method < 2/*vmsa::Method_stop*/; ++i_method)
          {
            string KEY_spectra_QA = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),vmsa::Dca_start,0/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            cout << KEY_spectra_QA << endl;
            h_mSpectra[KEY_spectra_QA];
            h_mSpectra[KEY_spectra_QA]->SetTitle("");
            h_mSpectra[KEY_spectra_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mSpectra[KEY_spectra_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mSpectra[KEY_spectra_QA]->GetXaxis()->SetTitle("p_{T}");
            h_mSpectra[KEY_spectra_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mSpectra[KEY_spectra_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mSpectra[KEY_spectra_QA]->GetXaxis()->CenterTitle();
  
            //h_mSpectra[KEY_spectra_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_spectra_QA]->GetMinimum(),1.2*h_mCounts[KEY_spectra_QA]->GetMaximum());
            h_mSpectra[KEY_spectra_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mSpectra[KEY_spectra_QA]->GetYaxis()->SetTitle("y");
            h_mSpectra[KEY_spectra_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mSpectra[KEY_spectra_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mSpectra[KEY_spectra_QA]->GetYaxis()->CenterTitle();
  
            h_mSpectra[KEY_spectra_QA]->Draw("colz");
          }
        }
        c_diff->SaveAs(outputname.c_str());
      }
    }
  }
#endif
  cout << "Done Plotting" << endl;

  //string outputfile = Form("../output/AuAu%s/%s/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string outputfile = Form("../output/AuAu%s/%s/TPCOnly_ptyspectra_RawPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/RawPhiPtSys_%s_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  //h_frame->Write();
  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
              for(int i_poly = 0; i_poly < 3/*3*/; i_poly++)
              {
	        string KEY_spectra = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	        h_mSpectra[KEY_spectra]->Write();
              }
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}
