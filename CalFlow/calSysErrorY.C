#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

//std::string centrality[13] = {"7080","6070","5060","4050","3040","2030","1020","0510","0005","0080","0010","1040","4080"};
//std::string centralityP[13] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%","0-80%","0-10%","10-40%","40-80%"};
std::string centrality[3] = {"4080","1040","0010"};
std::string centralityP[3] = {"40-80%","10-40%","0-10%"};

void calSysErrorY(Int_t energy = 4, Int_t pid = 0, int etamode = 0)
{
  gStyle->SetEndErrorSize(6);

  std::string etastring = "";
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";
  
  int ypadding = 0;
  //if(etamode == 0) ypadding = 2;
  //if(etamode == 3) ypadding = 5;
  //if(etamode == 4) ypadding = 4;
  //if(etamode == 5) ypadding = 3;

  string inputfile = Form("../output/AuAu%s/Flow/%s/RawPhiYSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mV2;
  TH1DMap h_mCounts;
  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++) // Centrality loop
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
                string KEY_v2 = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
                g_mV2[KEY_v2] = (TGraphAsymmErrors*)File_InPut->Get(KEY_v2.c_str());
                for(int i_y = 0+ypadding; i_y < vmsa::y_total-ypadding; i_y++) // Rapidity bin loop
                {
                  string KEY_counts = Form("y_%d_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_y,i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
                  h_mCounts[KEY_counts] = (TH1D*)File_InPut->Get(KEY_counts.c_str())->Clone();
                }
              }
            }
          }
        }
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPut->Get("h_frame");

  cout << "Loaded all the plots" << endl;
  
#if _PlotQA_
  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd();
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  //PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++) // Centrality loop
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
                string KEY_v2 = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
                Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mV2[KEY_v2],24,i_sigma+10*i_method+1,0,1.1);
              }
            }
          }
        }
      }
    }
  }
#endif
  TGraMap g_mSysErrors; 
  TGraMap g_mStatErrors; 
 
  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++) // Centrality loop
    {
      string KEY_Default = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

      g_mSysErrors[KEY_Default] = new TGraphAsymmErrors();
      g_mStatErrors[KEY_Default] = new TGraphAsymmErrors();
        
      double sysErr[g_mV2[KEY_Default]->GetN()][5]; // N points with 5 sources of systematics for each

      for(Int_t i_point = 0+ypadding; i_point < vmsa::y_total-ypadding; ++i_point)
      {
        double sysDca[9];
        double sysNSig[9];
        double sysNorm[9];
        //double sysSig[vmsa::Sig_stop];
        //double sysMeth[vmsa::Method_stop];

        double y_def, v2_def;
        g_mV2[KEY_Default]->GetPoint(i_point-ypadding,y_def,v2_def); 

        int idx = 0;
        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
        { 
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              if(i_dca == 0 && (i_sigma != 0 || i_method == 0)) continue;
              if(i_dca != 0 && i_sigma != 0 && i_method == 1) continue;

              string KEY = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
              double y_sys, v2_sys;
              g_mV2[KEY]->GetPoint(i_point-ypadding,y_sys,v2_sys);
              sysDca[idx] = v2_sys;
              idx++;
            }
          }
        }
        idx = 0;

        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              if(i_sig == 0 && (i_sigma != 0 || i_method == 0)) continue;
              if(i_sig != 0 && i_sigma != 0 && i_method == 1) continue;
              string KEY = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
              double y_sys, v2_sys;
              g_mV2[KEY]->GetPoint(i_point-ypadding,y_sys,v2_sys);
              sysNSig[idx] = v2_sys;
              idx++;
            }
          }
        }
        idx = 0;

        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              if(i_norm == 0 && (i_sigma != 0 || i_method == 0)) continue;
              if(i_norm != 0 && i_sigma != 0 && i_method == 1) continue;
              string KEY = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              double y_sys, v2_sys;
              g_mV2[KEY]->GetPoint(i_point-ypadding,y_sys,v2_sys);
              sysNorm[idx] = v2_sys;
              idx++;
            }
          }
        }
        idx = 0;	 

        /*for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
        {
          string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str());
          double y_sys, v2_sys;
          g_mV2[KEY]->GetPoint(i_point,y_sys,v2_sys);
          sysSig[i_sigma] = v2_sys;
        }*/

        /*for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str());
          double y_sys, v2_sys;
          g_mV2[KEY]->GetPoint(i_point,y_sys,v2_sys);
          sysMeth[i_method] = v2_sys;
        }*/

        Double_t v2_min[3] =  { TMath::MinElement(9,sysDca),
                                TMath::MinElement(9,sysNSig),
                                TMath::MinElement(9,sysNorm),     };

        Double_t v2_max[3] =  { TMath::MaxElement(9,sysDca),
                                TMath::MaxElement(9,sysNSig),
                                TMath::MaxElement(9,sysNorm),     };
      
        double SysError_v2 = 0.0;
        for(int i = 0; i < 3; i++)
        {
          double sourcei = TMath::Power((v2_max[i] - v2_min[i])/TMath::Sqrt(12.0),2);
          //cout << "v2_min = " << v2_min[i] << ", v2_max = " << v2_max[i] << endl;
          SysError_v2 += sourcei;
        }
 
        SysError_v2 = TMath::Sqrt(SysError_v2);

        Double_t pt, v2;
        g_mV2[KEY_Default]->GetPoint(i_point-ypadding,pt,v2);

        //float mean_v2 = total_v2/(float)counter;
        g_mSysErrors[KEY_Default]->SetPoint(i_point-ypadding,pt,v2);
        g_mSysErrors[KEY_Default]->SetPointError(i_point-ypadding,0.0,0.0,SysError_v2,SysError_v2);

        double StatError_v2 = g_mV2[KEY_Default]->GetErrorYhigh(i_point-ypadding);
        g_mStatErrors[KEY_Default]->SetPoint(i_point-ypadding,pt,v2);
        g_mStatErrors[KEY_Default]->SetPointError(i_point-ypadding,0.0,0.0,StatError_v2,StatError_v2);
      }
    }
  }

  cout << "Calculated the systematic errors" << endl;

  TGraphAsymmErrors *g_mSysErrorsInt = new TGraphAsymmErrors(); 
  TGraphAsymmErrors *g_mStatErrorsInt = new TGraphAsymmErrors(); 

  double weight_y = 0.0;
  double weight_error_stat_y = 0.0;
  double weight_error_sys_y = 0.0;
  double weight_all_y = 0.0;

  //for(int i_y = 0; i_y < vmsa::eta_total; i_y++) // Rapidity bin loop
  for(int i_y = 0+ypadding; i_y < vmsa::y_total-ypadding; i_y++) // Rapidity bin loop
  {
    //string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
    double weight_rho00 = 0;
    double weight_error_stat_rho00 = 0;
    double weight_error_sys_rho00 = 0;
    double weight_all = 0;

    for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
    {           
      //for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
      for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++) // Centrality loop
      { 
        string KEY_Default = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
        string KEY_counts = Form("y_%d_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_y,i_pt,i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
        //h_mCounts[KEY_counts];
        TH1F *PtCos = (TH1F*)h_mCounts[KEY_counts]->Clone();
        double weight = PtCos->Integral(1,7); 
         
        double etaStat, rhoStat, rhoErrStat;
        g_mStatErrors[KEY_Default]->GetPoint(i_y-ypadding,etaStat,rhoStat);   
        if(rhoStat == -999.0) { cout << "skip this rapidity" << endl; continue;}
        rhoErrStat = g_mStatErrors[KEY_Default]->GetErrorYhigh(i_y-ypadding);   
        double etaSys, rhoSys, rhoErrSys;
        g_mSysErrors[KEY_Default]->GetPoint(i_y-ypadding,etaSys,rhoSys);   
        rhoErrSys = g_mSysErrors[KEY_Default]->GetErrorYhigh(i_y-ypadding);  

        weight_rho00 += rhoStat*weight;
        cout << "y = " << i_y << ", Cent = " << i_cent << ", pT = " << i_pt << ",    rho00 = " << rhoStat << ", weight = " << weight << endl;
        weight_error_stat_rho00 += rhoErrStat*rhoErrStat*weight*weight;
        cout << "y = " << i_y << ", Cent = " << i_cent << ", pT = " << i_pt << ",    rho00 stat error = " << rhoErrStat << ", weight = " << weight << endl;
        weight_error_sys_rho00 += rhoErrSys*rhoErrSys*weight*weight;
        cout << "y = " << i_y << ", Cent = " << i_cent << ", pT = " << i_pt << ",    rho00 sys error  = " << rhoErrSys  << ", weight = " << weight << endl;
        weight_all += weight;
 
        weight_y += rhoStat*weight;
        weight_error_stat_y += rhoErrStat*rhoErrStat*weight*weight;
        weight_error_sys_y += rhoErrSys*rhoErrSys*weight*weight;
        weight_all_y += weight;       

      } 
    }
    weight_rho00 /= weight_all;
    cout << "y = " << i_y << ",    weight_rho00 = " << weight_rho00 << endl;
    weight_error_stat_rho00 = TMath::Sqrt(weight_error_stat_rho00)/weight_all;
    weight_error_sys_rho00 = TMath::Sqrt(weight_error_sys_rho00)/weight_all;

    double y_mean = (vmsa::ystart[i_y] + vmsa::ystop[i_y])/2.0;
    g_mStatErrorsInt->SetPoint(i_y-ypadding, y_mean, weight_rho00);
    g_mStatErrorsInt->SetPointError(i_y-ypadding, 0.0, 0.0, weight_error_stat_rho00, weight_error_stat_rho00);
    g_mSysErrorsInt->SetPoint(i_y-ypadding, y_mean, weight_rho00);
    g_mSysErrorsInt->SetPointError(i_y-ypadding, 0.0, 0.0, weight_error_sys_rho00, weight_error_sys_rho00);  
  }
  weight_y /= weight_all_y;
  weight_error_stat_y = TMath::Sqrt(weight_error_stat_y)/weight_all_y;
  weight_error_sys_y = TMath::Sqrt(weight_error_sys_y)/weight_all_y;

  cout << "Calculated integrated values over centrality and pt" << endl;

  cout << "Total Integrated rho00 = " << weight_y << " +/- " << weight_error_stat_y << " stat. " << " +/- " << weight_error_sys_y << "sys." << endl; 

  //cout << "All good" << endl;

  TCanvas *c_v2_SysError = new TCanvas("c_v2_SysError","c_v2_SysError",600,10,800,800);
  c_v2_SysError->cd();
  c_v2_SysError->cd()->SetLeftMargin(0.15);
  c_v2_SysError->cd()->SetBottomMargin(0.15);
  c_v2_SysError->cd()->SetTicks(1,1);
  c_v2_SysError->cd()->SetGrid(0,0);
  h_frame->Draw("pE");
  
  int tableNum[4] = {337,43,141,239};
  /*if(i_cent > 8)
  {
    TFile *besi_cd = TFile::Open("../data/BESI/HEPData-ins1395151-v2-root.root");
    TDirectory *dir = (TDirectory*) besi_cd->Get(Form("Table %d;1",tableNum[i_cent-9]));
    dir->cd();
    TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1;1");  
    besi19->SetLineColor(kBlack);
    besi19->SetMarkerStyle(20);
    besi19->SetMarkerSize(1.3);
    besi19->SetMarkerColor(kBlack);
    besi19->SetLineColor(kBlack);
    besi19->Draw("pE same");

    string leg_besi = "#phi BES-I";
    Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,20,kBlack,1.3);
    plotTopLegend((char*)leg_besi.c_str(),0.65,0.24,0.03,1,0.0,42,0);
  
    Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,20,kRed,1.3);
    string leg_count = "#phi BES-II";
    plotTopLegend((char*)leg_count.c_str(),0.65,0.26,0.03,1,0.0,42,0);
  }*/
 
 
  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),0.0,0.275,0.04,1,0.0,42,0);
  //plotTopLegend((char*)centralityP[i_cent].c_str(),3.1,0.245,0.04,1,0.0,42,0);
  plotTopLegend((char*)"0-80%",0.0,0.245,0.04,1,0.0,42,0);

  g_mStatErrorsInt->SetMarkerStyle(20);
  g_mStatErrorsInt->SetMarkerColor(kRed);
  g_mStatErrorsInt->SetLineColor(kRed);
  g_mStatErrorsInt->SetMarkerSize(1.3);

  g_mSysErrorsInt->SetMarkerStyle(20);
  g_mSysErrorsInt->SetMarkerSize(1.3);
  g_mSysErrorsInt->SetMarkerColor(kRed);
  g_mSysErrorsInt->SetLineColor(kRed);
  //mg_SysErrors->Draw("pE || same");
  g_mStatErrorsInt->Draw("pE same");

  plotSysErrorsBox(g_mSysErrorsInt,kRed);
  //PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  c_v2_SysError->SaveAs(Form("figures/AuAu%s/%s/V2_Rapidity_SysErrors_%s.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str()));

  string OutPutFile = Form("../output/AuAu%s/Flow/%s/V2_Rapidity_SysErrors_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {           
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++)
    {
      string KEY_Default = Form("v2Raw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      TString StatErrorV2 = Form("g_v2_%s_%s_pt%d_cent%d_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),i_pt,i_cent);
      g_mStatErrors[KEY_Default]->SetName(StatErrorV2.Data());
      g_mStatErrors[KEY_Default]->Write();
      TString SysErrorV2 = Form("g_v2_%s_%s_pt%d_cent%d_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),i_pt,i_cent);
      g_mSysErrors[KEY_Default]->SetName(SysErrorV2.Data());
      g_mSysErrors[KEY_Default]->Write();
    }
  }
  File_OutPut->Close();
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
