#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig6f()
{
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rho_1st_stat = (TGraphAsymmErrors*)File_Input->Get("rho00_1stEP_pt_stat_200");
  TGraphAsymmErrors *g_rho_1st_sys  = (TGraphAsymmErrors*)File_Input->Get("rho00_1stEP_pt_sys_200");
  TGraphAsymmErrors *g_rho_2nd_stat = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_pt_stat_200");
  TGraphAsymmErrors *g_rho_2nd_sys  = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_pt_sys_200");

  string outputfile_1st = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig6f_1st.txt";
  ofstream file_output_1st;
  file_output_1st.open(outputfile_1st.c_str());

  file_output_1st << "$p_{T} (GeV/c)$" << endl;
  file_output_1st << "1" << endl;
  file_output_1st << "$\\rho_{00}$" << endl;
  file_output_1st << "$1^{st}$-order EP" << endl;
  file_output_1st << "no" << endl;
  file_output_1st << "none" << endl;
  file_output_1st << "none" << endl;
  file_output_1st << "2" << endl;
  file_output_1st << "stat." << endl;
  file_output_1st << "symmetric" << endl;
  file_output_1st << "syst." << endl;
  file_output_1st << "symmetric" << endl;

  // 1st
  for(int iPt = 0; iPt < g_rho_1st_stat->GetN(); ++iPt)
  {
    double pt, rho;
    g_rho_1st_stat->GetPoint(iPt,pt,rho);
    double errStat = g_rho_1st_stat->GetErrorYhigh(iPt);
    double errSyst = g_rho_1st_sys->GetErrorYhigh(iPt);

    if(errStat > 0)
    {
      cout << pt << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output_1st << pt << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output_1st << "***" << endl;
  file_output_1st.close();

  string outputfile_2nd = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig6f_2nd.txt";
  ofstream file_output_2nd;
  file_output_2nd.open(outputfile_2nd.c_str());

  file_output_2nd << "$p_{T} (GeV/c)$" << endl;
  file_output_2nd << "1" << endl;
  file_output_2nd << "$\\rho_{00}$" << endl;
  file_output_2nd << "$1^{st}$-order EP" << endl;
  file_output_2nd << "no" << endl;
  file_output_2nd << "none" << endl;
  file_output_2nd << "none" << endl;
  file_output_2nd << "2" << endl;
  file_output_2nd << "stat." << endl;
  file_output_2nd << "symmetric" << endl;
  file_output_2nd << "syst." << endl;
  file_output_2nd << "symmetric" << endl;

  // 2nd
  for(int iPt = 0; iPt < g_rho_2nd_stat->GetN(); ++iPt)
  {
    double pt, rho;
    g_rho_2nd_stat->GetPoint(iPt,pt,rho);
    double errStat = g_rho_2nd_stat->GetErrorYhigh(iPt);
    double errSyst = g_rho_2nd_sys->GetErrorYhigh(iPt);

    if(errStat > 0)
    {
      cout << pt << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output_2nd << pt << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output_2nd << "***" << endl;
  file_output_2nd.close();

}
