#include <TGraphAsymmErrors.h>
#include "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_pT_July01_2020.h"

void makeData_rho00PtKstar()
{
  // 11 GeV
  TGraphAsymmErrors *g_rhoPt_stat_11GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_11GeV->SetName("rho00_2ndEP_pt_stat_11");
  TGraphAsymmErrors *g_rhoPt_sys_11GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_11GeV->SetName("rho00_2ndEP_pt_sys_11");
  for(int i_pt = 0; i_pt < kstar_Npt_11; ++i_pt)
  {
    g_rhoPt_stat_11GeV->SetPoint(i_pt,kstar_pt_11[i_pt],kstar_rho00_pT_11[i_pt]);
    g_rhoPt_stat_11GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_11_stat[i_pt],kstar_rho00_pT_11_stat[i_pt]);
    g_rhoPt_sys_11GeV->SetPoint(i_pt,kstar_pt_11[i_pt],kstar_rho00_pT_11[i_pt]);
    g_rhoPt_sys_11GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_11_syst[i_pt],kstar_rho00_pT_11_syst[i_pt]);
  }

  // 14 GeV
  TGraphAsymmErrors *g_rhoPt_stat_14GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_14GeV->SetName("rho00_2ndEP_pt_stat_14");
  TGraphAsymmErrors *g_rhoPt_sys_14GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_14GeV->SetName("rho00_2ndEP_pt_sys_14");
  for(int i_pt = 0; i_pt < kstar_Npt_14; ++i_pt)
  {
    g_rhoPt_stat_14GeV->SetPoint(i_pt,kstar_pt_14[i_pt],kstar_rho00_pT_14[i_pt]);
    g_rhoPt_stat_14GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_14_stat[i_pt],kstar_rho00_pT_14_stat[i_pt]);
    g_rhoPt_sys_14GeV->SetPoint(i_pt,kstar_pt_14[i_pt],kstar_rho00_pT_14[i_pt]);
    g_rhoPt_sys_14GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_14_syst[i_pt],kstar_rho00_pT_14_syst[i_pt]);
  }

  // 19 GeV
  TGraphAsymmErrors *g_rhoPt_stat_19GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_19GeV->SetName("rho00_2ndEP_pt_stat_19");
  TGraphAsymmErrors *g_rhoPt_sys_19GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_19GeV->SetName("rho00_2ndEP_pt_sys_19");
  for(int i_pt = 0; i_pt < kstar_Npt_19; ++i_pt)
  {
    g_rhoPt_stat_19GeV->SetPoint(i_pt,kstar_pt_19[i_pt],kstar_rho00_pT_19[i_pt]);
    g_rhoPt_stat_19GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_19_stat[i_pt],kstar_rho00_pT_19_stat[i_pt]);
    g_rhoPt_sys_19GeV->SetPoint(i_pt,kstar_pt_19[i_pt],kstar_rho00_pT_19[i_pt]);
    g_rhoPt_sys_19GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_19_syst[i_pt],kstar_rho00_pT_19_syst[i_pt]);
  }

  // 27 GeV
  TGraphAsymmErrors *g_rhoPt_stat_27GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_27GeV->SetName("rho00_2ndEP_pt_stat_27");
  TGraphAsymmErrors *g_rhoPt_sys_27GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_27GeV->SetName("rho00_2ndEP_pt_sys_27");
  for(int i_pt = 0; i_pt < kstar_Npt_27; ++i_pt)
  {
    g_rhoPt_stat_27GeV->SetPoint(i_pt,kstar_pt_27[i_pt],kstar_rho00_pT_27[i_pt]);
    g_rhoPt_stat_27GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_27_stat[i_pt],kstar_rho00_pT_27_stat[i_pt]);
    g_rhoPt_sys_27GeV->SetPoint(i_pt,kstar_pt_27[i_pt],kstar_rho00_pT_27[i_pt]);
    g_rhoPt_sys_27GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_27_syst[i_pt],kstar_rho00_pT_27_syst[i_pt]);
  }

  // 39 GeV
  TGraphAsymmErrors *g_rhoPt_stat_39GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_39GeV->SetName("rho00_2ndEP_pt_stat_39");
  TGraphAsymmErrors *g_rhoPt_sys_39GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_39GeV->SetName("rho00_2ndEP_pt_sys_39");
  for(int i_pt = 0; i_pt < kstar_Npt_39; ++i_pt)
  {
    g_rhoPt_stat_39GeV->SetPoint(i_pt,kstar_pt_39[i_pt],kstar_rho00_pT_39[i_pt]);
    g_rhoPt_stat_39GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_39_stat[i_pt],kstar_rho00_pT_39_stat[i_pt]);
    g_rhoPt_sys_39GeV->SetPoint(i_pt,kstar_pt_39[i_pt],kstar_rho00_pT_39[i_pt]);
    g_rhoPt_sys_39GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_39_syst[i_pt],kstar_rho00_pT_39_syst[i_pt]);
  }

  // 54 GeV
  TGraphAsymmErrors *g_rhoPt_stat_54GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_54GeV->SetName("rho00_2ndEP_pt_stat_54");
  TGraphAsymmErrors *g_rhoPt_sys_54GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_54GeV->SetName("rho00_2ndEP_pt_sys_54");
  for(int i_pt = 0; i_pt < kstar_Npt_54; ++i_pt)
  {
    g_rhoPt_stat_54GeV->SetPoint(i_pt,kstar_pt_54[i_pt],kstar_rho00_pT_54[i_pt]);
    g_rhoPt_stat_54GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_54_stat[i_pt],kstar_rho00_pT_54_stat[i_pt]);
    g_rhoPt_sys_54GeV->SetPoint(i_pt,kstar_pt_54[i_pt],kstar_rho00_pT_54[i_pt]);
    g_rhoPt_sys_54GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_54_syst[i_pt],kstar_rho00_pT_54_syst[i_pt]);
  }

  // 200 GeV
  TGraphAsymmErrors *g_rhoPt_stat_200GeV = new TGraphAsymmErrors();
  g_rhoPt_stat_200GeV->SetName("rho00_2ndEP_pt_stat_200");
  TGraphAsymmErrors *g_rhoPt_sys_200GeV = new TGraphAsymmErrors();
  g_rhoPt_sys_200GeV->SetName("rho00_2ndEP_pt_sys_200");
  for(int i_pt = 0; i_pt < kstar_Npt_200; ++i_pt)
  {
    g_rhoPt_stat_200GeV->SetPoint(i_pt,kstar_pt_200[i_pt],kstar_rho00_pT_200[i_pt]);
    g_rhoPt_stat_200GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_200_stat[i_pt],kstar_rho00_pT_200_stat[i_pt]);
    g_rhoPt_sys_200GeV->SetPoint(i_pt,kstar_pt_200[i_pt],kstar_rho00_pT_200[i_pt]);
    g_rhoPt_sys_200GeV->SetPointError(i_pt,0.0,0.0,kstar_rho00_pT_200_syst[i_pt],kstar_rho00_pT_200_syst[i_pt]);
  }

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_pT_July01_2020.root","RECREATE");
  File_OutPut->cd();
  g_rhoPt_stat_11GeV->Write();
  g_rhoPt_sys_11GeV->Write();
  g_rhoPt_stat_14GeV->Write();
  g_rhoPt_sys_14GeV->Write();
  g_rhoPt_stat_19GeV->Write();
  g_rhoPt_sys_19GeV->Write();
  g_rhoPt_stat_27GeV->Write();
  g_rhoPt_sys_27GeV->Write();
  g_rhoPt_stat_39GeV->Write();
  g_rhoPt_sys_39GeV->Write();
  g_rhoPt_stat_54GeV->Write();
  g_rhoPt_sys_54GeV->Write();
  g_rhoPt_stat_200GeV->Write();
  g_rhoPt_sys_200GeV->Write();
  File_OutPut->Close();
}
