#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig9_kstar()
{
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_0020.root");
  // K* STAR
  //beam-energy dependence of kstar rho00 from STAR, pT: 1.0 - 1.5 GeV/c, 0-20%
  TGraphAsymmErrors *g_rhoKstar_stat       = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoCentKstar_0020_stat");
  TGraphAsymmErrors *g_rhoKstar_sys        = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoCentKstar_0020_sys");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig9_kstar.txt";
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
    double errSyst = g_rhoKstar_sys->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}
