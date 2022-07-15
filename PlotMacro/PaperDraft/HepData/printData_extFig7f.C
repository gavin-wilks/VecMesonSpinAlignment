#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig7f()
{
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_pT_May21_2021.root");
  TGraphAsymmErrors *g_rho_2nd_stat = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_pt_stat_54");
  TGraphAsymmErrors *g_rho_2nd_sys  = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_pt_sys_54");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig7f.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$p_{T} (GeV/c)$" << endl;
  file_output << "1" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "$1^{st}$-order EP" << endl;
  file_output << "no" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "2" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;
  file_output << "syst." << endl;
  file_output << "symmetric" << endl;

  for(int iPt = 0; iPt < g_rho_2nd_stat->GetN(); ++iPt)
  {
    double pt, rho;
    g_rho_2nd_stat->GetPoint(iPt,pt,rho);
    double errStat = g_rho_2nd_stat->GetErrorYhigh(iPt);
    double errSyst = g_rho_2nd_sys->GetErrorYhigh(iPt);

    if(errStat > 0)
    {
      cout << pt << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << pt << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;
  file_output.close();
}
