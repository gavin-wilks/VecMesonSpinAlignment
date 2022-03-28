#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_fig3_kstar()
{
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_sNN_June9_2021.root");
  TGraphAsymmErrors *g_rhoKstar_stat       = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_stat");
  TGraphAsymmErrors *g_rhoKstar_sys        = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_sys");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataFig3_kstar.txt";
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
  for(int iEnergy = 0; iEnergy < g_rhoKstar_stat->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoKstar_stat->GetPoint(iEnergy,energy,rho);
    double errStat = g_rhoKstar_stat->GetErrorYhigh(iEnergy);
    double errSyst = g_rhoKstar_stat->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}
