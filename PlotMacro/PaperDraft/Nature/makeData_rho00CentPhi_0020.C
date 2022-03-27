#include <TGraphAsymmErrors.h>
#include "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rhoCent_0020_LXaxis.h"

void makeData_rho00CentPhi_0020()
{
  TGraphAsymmErrors *g_rho1st_0020_stat = new TGraphAsymmErrors();
  g_rho1st_0020_stat->SetName("g_rho1st_0020_stat");
  TGraphAsymmErrors *g_rho1st_0020_sys = new TGraphAsymmErrors();
  g_rho1st_0020_sys->SetName("g_rho1st_0020_sys");

  TGraphAsymmErrors *g_rho2nd_0020_stat = new TGraphAsymmErrors();
  g_rho2nd_0020_stat->SetName("g_rho2nd_0020_stat");
  TGraphAsymmErrors *g_rho2nd_0020_sys = new TGraphAsymmErrors();
  g_rho2nd_0020_sys->SetName("g_rho2nd_0020_sys");

  for(int iEnergy = 0; iEnergy < 4; ++iEnergy)
  {
    g_rho1st_0020_stat->SetPoint(iEnergy,energy[iEnergy]*1.08,rho_1st[iEnergy]);
    g_rho1st_0020_stat->SetPointError(iEnergy,0.0,0.0,errStat_1st[iEnergy],errStat_1st[iEnergy]);
    g_rho1st_0020_sys->SetPoint(iEnergy,energy[iEnergy]*1.08,rho_1st[iEnergy]);
    g_rho1st_0020_sys->SetPointError(iEnergy,0.0,0.0,errSys_1st[iEnergy],errSys_1st[iEnergy]);

    g_rho2nd_0020_stat->SetPoint(iEnergy,energy[iEnergy],rho_2nd[iEnergy]);
    g_rho2nd_0020_stat->SetPointError(iEnergy,0.0,0.0,errStat_2nd[iEnergy],errStat_2nd[iEnergy]);
    g_rho2nd_0020_sys->SetPoint(iEnergy,energy[iEnergy],rho_2nd[iEnergy]);
    g_rho2nd_0020_sys->SetPointError(iEnergy,0.0,0.0,errSys_2nd[iEnergy],errSys_2nd[iEnergy]);
  }


  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rhoCent_0020_LXaxis.root","RECREATE");
  File_OutPut->cd();
  g_rho1st_0020_stat->Write();
  g_rho1st_0020_sys->Write();

  g_rho2nd_0020_stat->Write();
  g_rho2nd_0020_sys->Write();
  File_OutPut->Close();
}
