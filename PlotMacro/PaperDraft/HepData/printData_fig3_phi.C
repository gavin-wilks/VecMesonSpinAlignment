#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_fig3_phi()
{
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_sys");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataFig3_phi.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$\\sqrt{s_{NN}} (GeV)$" << endl;
  file_output << "1" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "2nd EP" << endl;
  file_output << "no" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "2" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;
  file_output << "syst." << endl;
  file_output << "symmetric" << endl;

  // 2nd EP
  for(int iEnergy = 0; iEnergy < g_rhoPhi_2nd_stat_Laxis->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Laxis->GetPoint(iEnergy,energy,rho);
    double errStat = g_rhoPhi_2nd_stat_Laxis->GetErrorYhigh(iEnergy);
    double errSyst = g_rhoPhi_2nd_sys_Laxis->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}
