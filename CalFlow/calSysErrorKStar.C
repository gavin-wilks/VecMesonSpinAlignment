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

std::string centrality[13] = {"7080","6070","5060","4050","3040","2030","1020","0510","0005","0080","0010","1040","4080"};
std::string centralityP[13] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%","0-80%","0-10%","10-40%","40-80%"};

void calSysErrorKStar(Int_t energy = 4, Int_t pid = 2, Int_t i_cent = 9)
{
  gStyle->SetEndErrorSize(6);
  
  string inputfile = Form("../output/AuAu%s/Flow/%s/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mV2;
  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      if( i_dca != 0 && i_sig != 0 ) continue;
      for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
      {
        if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              string KEY_v2 = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              g_mV2[KEY_v2] = (TGraphAsymmErrors*)File_InPut->Get(KEY_v2.c_str());
            }
          }
        }
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPut->Get("h_frame");

#if _PlotQA_
  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd();
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  //PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      if( i_dca != 0 && i_sig != 0 ) continue;
      for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
      {
        if(( i_dca != 0 && i_nhit != 0 ) || ( i_sig != 0 && i_nhit != 0 )) continue;
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              string KEY_v2 = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mV2[KEY_v2],24,i_sigma+10*i_method+1,1.1);
            }
          }
        }
      }
    }
  }
#endif
 
  string KEY_Default = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();
    
  double sysErr[g_mV2[KEY_Default]->GetN()][6]; // N points with 6 sources of systematics for each

  for(Int_t i_point = 0; i_point < g_mV2[KEY_Default]->GetN(); ++i_point)
  {
    double sysDca[vmsa::Dca_stop];
    double sysNSig[vmsa::nSigKaon_stop];
    double sysNHit[vmsa::mNHit_stop];
    double sysNorm[vmsa::Norm_stop];
    double sysSig[vmsa::Sig_stop];
    double sysMeth[vmsa::Method_stop];

    double pt_def, v2_def;
    g_mV2[KEY_Default]->GetPoint(i_point,pt_def,v2_def); 

    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    { 
      string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, v2_sys;
      g_mV2[KEY]->GetPoint(i_point,pt_sys,v2_sys);
      sysDca[i_dca] = v2_sys;
    }

    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,i_sig,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, v2_sys;
      g_mV2[KEY]->GetPoint(i_point,pt_sys,v2_sys);
      sysNSig[i_sig] = v2_sys;
    }

    for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
    {
      string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,i_nhit,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, v2_sys;
      g_mV2[KEY]->GetPoint(i_point,pt_sys,v2_sys);
      sysNHit[i_sig] = v2_sys;
    }

    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
    {
      string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,0,vmsa::mPID[pid].c_str(),i_norm,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, v2_sys;
      g_mV2[KEY]->GetPoint(i_point,pt_sys,v2_sys);
      sysNorm[i_norm] = v2_sys;
    }	 

    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str());
      double pt_sys, v2_sys;
      g_mV2[KEY]->GetPoint(i_point,pt_sys,v2_sys);
      sysSig[i_sigma] = v2_sys;
    }

    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
    {
      string KEY = Form("v2Raw_Centrality_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str());
      double pt_sys, v2_sys;
      g_mV2[KEY]->GetPoint(i_point,pt_sys,v2_sys);
      sysMeth[i_method] = v2_sys;
    }

    Double_t v2_min[6] = { TMath::MinElement(vmsa::Dca_stop,sysDca),
                           TMath::MinElement(vmsa::nSigKaon_stop,sysNSig),
                           TMath::MinElement(vmsa::mNHit_stop,sysNHit),
                           TMath::MinElement(vmsa::Norm_stop,sysNorm),
                           TMath::MinElement(vmsa::Sig_start,sysSig),
                           TMath::MinElement(vmsa::Method_stop,sysMeth)    };

    Double_t v2_max[6] = { TMath::MaxElement(vmsa::Dca_stop,sysDca),
                           TMath::MaxElement(vmsa::nSigKaon_stop,sysNSig),
                           TMath::MaxElement(vmsa::mNHit_stop,sysNHit),
                           TMath::MaxElement(vmsa::Norm_stop,sysNorm),
                           TMath::MaxElement(vmsa::Sig_start,sysSig),
                           TMath::MaxElement(vmsa::Method_stop,sysMeth)    };
  
    double SysError_v2 = 0.0;
    for(int i = 0; i < 6; i++)
    {
      double sourcei = TMath::Power((v2_max[i] - v2_min[i])/TMath::Sqrt(12.0),2);
      //cout << "v2_min = " << v2_min[i] << ", v2_max = " << v2_max[i] << endl;
      SysError_v2 += sourcei;
    }
 
    SysError_v2 = TMath::Sqrt(SysError_v2);

    Double_t pt, v2;
    g_mV2[KEY_Default]->GetPoint(i_point,pt,v2);

    //float mean_v2 = total_v2/(float)counter;
    g_SysErrors->SetPoint(i_point,pt,v2);
    g_SysErrors->SetPointError(i_point,0.0,0.0,SysError_v2,SysError_v2);

    double StatError_v2 = g_mV2[KEY_Default]->GetErrorYhigh(i_point);
    g_StatErrors->SetPoint(i_point,pt,v2);
    g_StatErrors->SetPointError(i_point,0.0,0.0,StatError_v2,StatError_v2);
  }

  //cout << "All good" << endl;

  TCanvas *c_v2_SysError = new TCanvas("c_v2_SysError","c_v2_SysError",600,10,800,800);
  c_v2_SysError->cd();
  c_v2_SysError->cd()->SetLeftMargin(0.15);
  c_v2_SysError->cd()->SetBottomMargin(0.15);
  c_v2_SysError->cd()->SetTicks(1,1);
  c_v2_SysError->cd()->SetGrid(0,0);
  h_frame->Draw("pE");
  
  int tableNum[4] = {337,43,141,239};
  if(i_cent > 8)
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
    besi19->Draw("pE Z same");

    string leg_besi = "#phi BES-I";
    Draw_TGAE_Point_new_Symbol(0.5,0.24,0.0,0.0,0.0,0.0,20,kBlack,1.3);
    plotTopLegend((char*)leg_besi.c_str(),0.65,0.24,0.03,1,0.0,42,0);
  
    Draw_TGAE_Point_new_Symbol(0.5,0.26,0.0,0.0,0.0,0.0,20,kRed,1.3);
    string leg_count = "#phi BES-II";
    plotTopLegend((char*)leg_count.c_str(),0.65,0.26,0.03,1,0.0,42,0);
  }
 
 
  string leg_energy = Form("AuAu %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),2.5,0.275,0.04,1,0.0,42,0);
  plotTopLegend((char*)centralityP[i_cent].c_str(),3.1,0.245,0.04,1,0.0,42,0);

  g_StatErrors->SetMarkerStyle(20);
  g_StatErrors->SetMarkerColor(kRed);
  g_StatErrors->SetLineColor(kRed);
  g_StatErrors->SetMarkerSize(1.3);

  g_SysErrors->SetMarkerStyle(20);
  g_SysErrors->SetMarkerSize(1.3);
  g_SysErrors->SetMarkerColor(kRed);
  g_SysErrors->SetLineColor(kRed);
  //g_SysErrors->Draw("pE || same");
  g_StatErrors->Draw("pE Z same");

  plotSysErrorsBox(g_SysErrors,kBlack,energy);
  //PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  c_v2_SysError->SaveAs(Form("figures/AuAu%s/%s/V2_SysErrors_%s.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),centrality[i_cent].c_str()));

  string OutPutFile = Form("../output/AuAu%s/Flow/%s/V2_SysErrors_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),centrality[i_cent].c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  TString StatErrorV2 = Form("g_v2_%s_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),centrality[i_cent].c_str());
  g_StatErrors->SetName(StatErrorV2.Data());
  g_StatErrors->Write();
  TString SysErrorV2 = Form("g_v2_%s_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),centrality[i_cent].c_str());
  g_SysErrors->SetName(SysErrorV2.Data());
  g_SysErrors->Write();
  File_OutPut->Close();
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
