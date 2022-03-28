#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig1c()
{
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/Second_ActualFit_27GeV_pT_2.root");
  TH1D *h_mMassPhi = (TH1D*)File_InputPhi->Get("imass_all_pt2")->Clone("h_mMassPhi");
  TH1D *h_mYieldsPhi = (TH1D*)File_InputPhi->Get("raw_yield")->Clone("h_mYieldsPhi");
  h_mYieldsPhi->Draw("pE");

  const long numOfEvent = 229186355;
  const double binWidth = h_mMassPhi->GetBinWidth(1);
  h_mYieldsPhi->Scale(1.0/(numOfEvent*binWidth));

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig1c.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$\\left|cos(\\theta*)\\right|$" << endl;
  file_output << "1" << endl;
  file_output << "Raw Yields" << endl;
  file_output << "$20-60\%$ centrality" << endl;
  file_output << "yes" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "1" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;

  for(int iBin = 0; iBin < h_mYieldsPhi->GetNbinsX(); ++iBin)
  {
    double binLowEdge  = h_mYieldsPhi->GetBinLowEdge(iBin+1);
    double binWitdh    = h_mYieldsPhi->GetBinWidth(iBin+1);
    double binHighEdge = binLowEdge+binWitdh;
    double binContent  = h_mYieldsPhi->GetBinContent(iBin+1);
    double binError    = h_mYieldsPhi->GetBinError(iBin+1);
    if(binError > 0)
    {
      // cout << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
      file_output << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}
