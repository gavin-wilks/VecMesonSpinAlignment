#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>

using namespace std;

void plotQA_TofMatchingEfficiency(int mEnergy = 6, int mPID = 0)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  int const mPlotStyle[10] = {0,0,1,1,1,1,0,0,0,1};

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch_2060.root",mBeamEnergy[mEnergy].c_str(),mBeamEnergy[mEnergy].c_str());
  cout << "open input file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1D *h_mEfficiency[10]; // 9: 20-60%, 0-8 from StRefMultCorr
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName = Form("h_mEfficiency_%s_Cent_%d",mParType[mPID].c_str(),i_cent);
    h_mEfficiency[i_cent] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
  }

  TH1D *h_play = new TH1D("h_play","h_play",81,-0.01,8.00);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,0.1);
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetGrid(0,0);
  c_play->cd()->SetTicks(1,1);

  string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  h_play->SetTitle(title.c_str());
  h_play->SetStats(0);

  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetTitle("Efficiency");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetRangeUser(0.0,1.05);
  h_play->Draw("h");

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(mPlotStyle[i_cent] > 0)
    {
      h_mEfficiency[i_cent]->SetLineColor(i_cent+1);
      if(i_cent == 9) h_mEfficiency[i_cent]->SetLineColor(2);
      h_mEfficiency[i_cent]->DrawCopy("HIST same");
    }
  }

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(mPlotStyle[i_cent] > 0)
    {
      leg->AddEntry(h_mEfficiency[i_cent],Centrality[i_cent].c_str());
    }
  }
  leg->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/TPC/c_TpcEffCentCom_%s%s.eps",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  // c_play->SaveAs(FigName.c_str());

  /*
  const int cent9 = 4;
  TH1D *h_mEfficiencyEta[12]; // phi-bin
  */
}
