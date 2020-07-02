#include <TGraphAsymmErrors.h>
#include "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_July01_2020.h"

void makeData_rho00CentKstar()
{
  // 39 GeV
  TGraphAsymmErrors *g_rhoCent_stat_39GeV = new TGraphAsymmErrors();
  g_rhoCent_stat_39GeV->SetName("g_rhoCent_39GeV_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_sys_39GeV = new TGraphAsymmErrors();
  g_rhoCent_sys_39GeV->SetName("g_rhoCent_39GeV_2nd_sys");
  for(int i_cent = 0; i_cent < kstar_Ncent39; ++i_cent)
  {
    g_rhoCent_stat_39GeV->SetPoint(i_cent,kstar_Centrality_39[i_cent],kstar_rho00_Cent_39[i_cent]);
    g_rhoCent_stat_39GeV->SetPointError(i_cent,0.0,0.0,kstar_rho00_Cent_39_stat[i_cent],kstar_rho00_Cent_39_stat[i_cent]);
    g_rhoCent_sys_39GeV->SetPoint(i_cent,kstar_Centrality_39[i_cent],kstar_rho00_Cent_39[i_cent]);
    g_rhoCent_sys_39GeV->SetPointError(i_cent,0.0,0.0,kstar_rho00_Cent_39_syst[i_cent],kstar_rho00_Cent_39_syst[i_cent]);
  }

  // 54 GeV
  TGraphAsymmErrors *g_rhoCent_stat_54GeV = new TGraphAsymmErrors();
  g_rhoCent_stat_54GeV->SetName("g_rhoCent_54GeV_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_sys_54GeV = new TGraphAsymmErrors();
  g_rhoCent_sys_54GeV->SetName("g_rhoCent_54GeV_2nd_sys");
  for(int i_cent = 0; i_cent < kstar_Ncent54; ++i_cent)
  {
    g_rhoCent_stat_54GeV->SetPoint(i_cent,kstar_Centrality_54[i_cent],kstar_rho00_Cent_54[i_cent]);
    g_rhoCent_stat_54GeV->SetPointError(i_cent,0.0,0.0,kstar_rho00_Cent_54_stat[i_cent],kstar_rho00_Cent_54_stat[i_cent]);
    g_rhoCent_sys_54GeV->SetPoint(i_cent,kstar_Centrality_54[i_cent],kstar_rho00_Cent_54[i_cent]);
    g_rhoCent_sys_54GeV->SetPointError(i_cent,0.0,0.0,kstar_rho00_Cent_54_syst[i_cent],kstar_rho00_Cent_54_syst[i_cent]);
  }

  // 200 GeV
  TGraphAsymmErrors *g_rhoCent_stat_200GeV = new TGraphAsymmErrors();
  g_rhoCent_stat_200GeV->SetName("g_rhoCent_200GeV_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_sys_200GeV = new TGraphAsymmErrors();
  g_rhoCent_sys_200GeV->SetName("g_rhoCent_200GeV_2nd_sys");
  for(int i_cent = 0; i_cent < kstar_Ncent200; ++i_cent)
  {
    g_rhoCent_stat_200GeV->SetPoint(i_cent,kstar_Centrality_200[i_cent],kstar_rho00_Cent_200[i_cent]);
    g_rhoCent_stat_200GeV->SetPointError(i_cent,0.0,0.0,kstar_rho00_Cent_200_stat[i_cent],kstar_rho00_Cent_200_stat[i_cent]);
    g_rhoCent_sys_200GeV->SetPoint(i_cent,kstar_Centrality_200[i_cent],kstar_rho00_Cent_200[i_cent]);
    g_rhoCent_sys_200GeV->SetPointError(i_cent,0.0,0.0,kstar_rho00_Cent_200_syst[i_cent],kstar_rho00_Cent_200_syst[i_cent]);
  }

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_July01_2020.root","RECREATE");
  File_OutPut->cd();
  g_rhoCent_stat_39GeV->Write();
  g_rhoCent_sys_39GeV->Write();
  g_rhoCent_stat_54GeV->Write();
  g_rhoCent_sys_54GeV->Write();
  g_rhoCent_stat_200GeV->Write();
  g_rhoCent_sys_200GeV->Write();
  File_OutPut->Close();
}
