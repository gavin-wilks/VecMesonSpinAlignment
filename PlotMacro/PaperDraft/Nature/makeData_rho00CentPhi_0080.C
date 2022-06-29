#include <TGraphAsymmErrors.h>
#include "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/NewF_JHChen/rhoCent_small_bin_Laxis.h"

void makeData_rho00CentPhi_0080()
{
  // 62GeV
  TGraphAsymmErrors *g_rhoCent_62GeV_1st_stat = new TGraphAsymmErrors();
  g_rhoCent_62GeV_1st_stat->SetName("g_rhoCent_62GeV_1st_stat");
  TGraphAsymmErrors *g_rhoCent_62GeV_1st_sys = new TGraphAsymmErrors();
  g_rhoCent_62GeV_1st_sys->SetName("g_rhoCent_62GeV_1st_sys");

  TGraphAsymmErrors *g_rhoCent_62GeV_2nd_stat = new TGraphAsymmErrors();
  g_rhoCent_62GeV_2nd_stat->SetName("g_rhoCent_62GeV_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_62GeV_2nd_sys = new TGraphAsymmErrors();
  g_rhoCent_62GeV_2nd_sys->SetName("g_rhoCent_62GeV_2nd_sys");

  for(int iCent = 0; iCent < 9; ++iCent)
  {
    g_rhoCent_62GeV_1st_stat->SetPoint(iCent,cent_1st[iCent],rho_1st_62[iCent]);
    g_rhoCent_62GeV_1st_stat->SetPointError(iCent,0.0,0.0,errStat_1st_62[iCent],errStat_1st_62[iCent]);
    g_rhoCent_62GeV_1st_sys->SetPoint(iCent,cent_1st[iCent],rho_1st_62[iCent]);
    g_rhoCent_62GeV_1st_sys->SetPointError(iCent,0.0,0.0,errSys_1st_62[iCent],errSys_1st_62[iCent]);

    g_rhoCent_62GeV_2nd_stat->SetPoint(iCent,cent_2nd[iCent],rho_2nd_62[iCent]);
    g_rhoCent_62GeV_2nd_stat->SetPointError(iCent,0.0,0.0,errStat_2nd_62[iCent],errStat_2nd_62[iCent]);
    g_rhoCent_62GeV_2nd_sys->SetPoint(iCent,cent_2nd[iCent],rho_2nd_62[iCent]);
    g_rhoCent_62GeV_2nd_sys->SetPointError(iCent,0.0,0.0,errSys_2nd_62[iCent],errSys_2nd_62[iCent]);
  }

  TFile *File_OutPut_62GeV = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/NewF_JHChen/rhoCent_62GeV_LXaxis.root","RECREATE");
  File_OutPut_62GeV->cd();
  g_rhoCent_62GeV_1st_stat->Write();
  g_rhoCent_62GeV_1st_sys->Write();

  g_rhoCent_62GeV_2nd_stat->Write();
  g_rhoCent_62GeV_2nd_sys->Write();
  File_OutPut_62GeV->Close();

  // 39GeV
  TGraphAsymmErrors *g_rhoCent_39GeV_1st_stat = new TGraphAsymmErrors();
  g_rhoCent_39GeV_1st_stat->SetName("g_rhoCent_39GeV_1st_stat");
  TGraphAsymmErrors *g_rhoCent_39GeV_1st_sys = new TGraphAsymmErrors();
  g_rhoCent_39GeV_1st_sys->SetName("g_rhoCent_39GeV_1st_sys");

  TGraphAsymmErrors *g_rhoCent_39GeV_2nd_stat = new TGraphAsymmErrors();
  g_rhoCent_39GeV_2nd_stat->SetName("g_rhoCent_39GeV_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_39GeV_2nd_sys = new TGraphAsymmErrors();
  g_rhoCent_39GeV_2nd_sys->SetName("g_rhoCent_39GeV_2nd_sys");

  for(int iCent = 0; iCent < 9; ++iCent)
  {
    g_rhoCent_39GeV_1st_stat->SetPoint(iCent,cent_1st[iCent],rho_1st_39[iCent]);
    g_rhoCent_39GeV_1st_stat->SetPointError(iCent,0.0,0.0,errStat_1st_39[iCent],errStat_1st_39[iCent]);
    g_rhoCent_39GeV_1st_sys->SetPoint(iCent,cent_1st[iCent],rho_1st_39[iCent]);
    g_rhoCent_39GeV_1st_sys->SetPointError(iCent,0.0,0.0,errSys_1st_39[iCent],errSys_1st_39[iCent]);

    g_rhoCent_39GeV_2nd_stat->SetPoint(iCent,cent_2nd[iCent],rho_2nd_39[iCent]);
    g_rhoCent_39GeV_2nd_stat->SetPointError(iCent,0.0,0.0,errStat_2nd_39[iCent],errStat_2nd_39[iCent]);
    g_rhoCent_39GeV_2nd_sys->SetPoint(iCent,cent_2nd[iCent],rho_2nd_39[iCent]);
    g_rhoCent_39GeV_2nd_sys->SetPointError(iCent,0.0,0.0,errSys_2nd_39[iCent],errSys_2nd_39[iCent]);
  }

  TFile *File_OutPut_39GeV = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/NewF_JHChen/rhoCent_39GeV_LXaxis.root","RECREATE");
  File_OutPut_39GeV->cd();
  g_rhoCent_39GeV_1st_stat->Write();
  g_rhoCent_39GeV_1st_sys->Write();

  g_rhoCent_39GeV_2nd_stat->Write();
  g_rhoCent_39GeV_2nd_sys->Write();
  File_OutPut_39GeV->Close();

  // 27GeV Run18
  TGraphAsymmErrors *g_rhoCent_27GeV_Run18_1st_stat = new TGraphAsymmErrors();
  g_rhoCent_27GeV_Run18_1st_stat->SetName("g_rhoCent_27GeV_Run18_1st_stat");
  TGraphAsymmErrors *g_rhoCent_27GeV_Run18_1st_sys = new TGraphAsymmErrors();
  g_rhoCent_27GeV_Run18_1st_sys->SetName("g_rhoCent_27GeV_Run18_1st_sys");

  TGraphAsymmErrors *g_rhoCent_27GeV_Run18_2nd_stat = new TGraphAsymmErrors();
  g_rhoCent_27GeV_Run18_2nd_stat->SetName("g_rhoCent_27GeV_Run18_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_27GeV_Run18_2nd_sys = new TGraphAsymmErrors();
  g_rhoCent_27GeV_Run18_2nd_sys->SetName("g_rhoCent_27GeV_Run18_2nd_sys");

  for(int iCent = 0; iCent < 9; ++iCent)
  {
    g_rhoCent_27GeV_Run18_1st_stat->SetPoint(iCent,cent_1st[iCent],rho_1st_27[iCent]);
    g_rhoCent_27GeV_Run18_1st_stat->SetPointError(iCent,0.0,0.0,errStat_1st_27[iCent],errStat_1st_27[iCent]);
    g_rhoCent_27GeV_Run18_1st_sys->SetPoint(iCent,cent_1st[iCent],rho_1st_27[iCent]);
    g_rhoCent_27GeV_Run18_1st_sys->SetPointError(iCent,0.0,0.0,errSys_1st_27[iCent],errSys_1st_27[iCent]);

    g_rhoCent_27GeV_Run18_2nd_stat->SetPoint(iCent,cent_2nd[iCent],rho_2nd_27[iCent]);
    g_rhoCent_27GeV_Run18_2nd_stat->SetPointError(iCent,0.0,0.0,errStat_2nd_27[iCent],errStat_2nd_27[iCent]);
    g_rhoCent_27GeV_Run18_2nd_sys->SetPoint(iCent,cent_2nd[iCent],rho_2nd_27[iCent]);
    g_rhoCent_27GeV_Run18_2nd_sys->SetPointError(iCent,0.0,0.0,errSys_2nd_27[iCent],errSys_2nd_27[iCent]);
  }

  TFile *File_OutPut_27GeV_Run18 = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/NewF_JHChen/rhoCent_27GeV_Run18_LXaxis.root","RECREATE");
  File_OutPut_27GeV_Run18->cd();
  g_rhoCent_27GeV_Run18_1st_stat->Write();
  g_rhoCent_27GeV_Run18_1st_sys->Write();

  g_rhoCent_27GeV_Run18_2nd_stat->Write();
  g_rhoCent_27GeV_Run18_2nd_sys->Write();
  File_OutPut_27GeV_Run18->Close();
}
